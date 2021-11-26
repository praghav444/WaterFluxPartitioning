#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Water Flux Partitioning using underlying Water Use Efficiency (uWUE) Method
# Writen by: Pushpendra Raghav on July 27, 2021
# @ppushpendra@crimson.ua.edu
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# uWUE <- GPP*sqr(VPD)/ET      -------(1)
# T:ET <- uWUEa/uWUEp
# where uWUEa is Apparent uWUE [estimated as the linear regression slope from a 
# moving window spanning either one or eight days (depending on desired 
# smoothening or data availability) or directly from Eqn 1 when estimating at 
# 30 minutes temporal resolution]
# uWUEp is potential uWUE [calculated at annual or seasonal scale using 95th 
# percentile regression between GPP*sqrt(VPD) and ET]
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Some important functions

# LE to ET
LE2ET <- function(LE, Ta){
  "
  Convert LE to ET
  
  Parameters
  ---------------
  LE: Latent Heat flux (W m-2)
  Ta: Air temperature (deg C)
  
  Returns
  --------------
  ET: Evapotranspiration (mm s-1)
  "
  lambda <- (2.501 - 0.00237*Ta)*1e6  # Latent heat of vaporization (J kg-1)
  ET <- LE/lambda
  return(ET)
}


calc_PET_PT <- function(Ta, Pa, Rn, G=NA, S=NA, alpha=NA){
  "
  Calculates potential evapotranspiration (PET) according to Priestley & Taylor 1972
  LE_pot = alpha * delta * (Rn - G)) / (delta + gamma)}
  
  Parameters
    ----------
    Ta : Air temperature (deg C)
    Pa : Atmospheric pressure (kPa)
    Rn : Net radiation (W m-2)
    G : Ground heat flux (W m-2); optional
    S : Sum of all storage fluxes (W m-2); optional
    alpha : Priestley-Taylor coefficient (default = 1.26)
    
  Returns
    -------
    PET : Potential evapotranspiration (kg m-2 s-1 or mm/s)
    LE_pot : Potential latent heat flux (W m-2)
    
  
  References
    ----------
    - Priestley, C.H.B., Taylor, R.J., 1972: On the assessment of surface heat flux
      and evaporation using large-scale parameters. Monthly Weather Review 100, 81-92.
  "
  G[is.na(G)] <- 0  # Set G to 0 if not provided
  S[is.na(S)] <- 0  # Set S to 0 if not provided
  lambda <- (2.501 - 0.00237*Ta)*1e6  # Latent heat of vaporization (J kg-1)
  Cp <- 1004.834   # specific heat of air for constant pressure (J K-1 kg-1)
  eps <- 0.622     # ratio of the molecular weight of water vapor to dry air (=Mw/Md)
  gamma  <- (Cp * Pa) / (eps * lambda)  # Psychrometric constant (kPa K-1)
  esat <- 0.6108*exp(17.27*Ta/(Ta+237.3))*1e3 # Saturation vapor pressure (Pa)
  delta <- 4098*esat/(237.3+Ta)^2*1e-3 # Slope of the saturation vapor pressure curve (kPa K-1)
  alpha[is.na(alpha)] <- 1.26
  LE_pot = (alpha * delta * (Rn - G - S)) / (delta + gamma)
  PET = LE2ET(LE_pot, Ta)
  return(data.frame(LE_pot, PET))
}

quantile_reg <- function(x, y , PolyDeg=1, tau=0.95, weights){
  "
  Quantile regression
    Fits a polynomial function (of degree PolyDeg) using quantile regression based on a percentile (tau).
    Based on script by Dr. Phillip M. Feldman, and based on method by Koenker, Roger, and
    Gilbert Bassett Jr. Regression Quantiles. Econometrica: Journal of the Econometric Society, 1978, 33-50.
  
  Parameters
    ----------
    x : independent variable
    y : dependent variable
    PolyDeg : Degree of polynomial function
    tau : Percentile for the data to fit to [0-1]
    weights : Vector to weight each point, must be same size as x
    
  Returns
    -------
    The resulting parameters in order of degree from high to low
  
  "
  model <- function(x, beta){
    "
    This example defines the model as a polynomial, where the coefficients of the
    polynomial are passed via `beta`.
    "
    library(signal)
    if(PolyDeg==0){
      return(x*beta)
    } else{
      return(polyval(beta, x))
    }
  }
  N_coefficients <- PolyDeg+1
  
  tilted_abs <- function(tau, x, weights){
    "
     The tilted absolute value function is used in quantile regression.
      INPUTS
       tau: This parameter is a probability, and thus takes values between 0 and 1.
       x: This parameter represents a value of the independent variable, and in
       general takes any real value (float).
    "
    return (weights * x * (tau - (x < 0)))
  }
  objective <- function(beta, tau, weights){
    "
    The objective function to be minimized is the sum of the tilted absolute
    values of the differences between the observations and the model.
    "
    return(sum(tilted_abs(tau, y - model(x, beta), weights)))
  }
  # Build weights if they don't exits:
  if(is.na(weights)){
    weights <- rep(1,length(x))
  }
  # Define starting point for optimization:
  beta_0 <- rep(0, N_coefficients)  
  if(N_coefficients >= 2){
    beta_0[1] <- 1.0
  }
  # `beta_hat[i]` will store the parameter estimates for the quantile
  # corresponding to `fractions[i]`:
  library(stats)
  beta_hat <- optim(beta_0, objective,weights=weights, tau=tau)
  beta_hat <- beta_hat$par
  return(beta_hat)
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Data Needed: GPP, VPD, ET
#-------------------------------------------------------------------------------
library(openxlsx)
library(lubridate)
library(dplyr)
library(data.table)
library(quantreg)
#-------------------------------------------------------------------------------
# Creating a dataframe having all the required data
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
df_met <- read.xlsx('C://Users/ppushpendra.STUDENT/Box/FVS_theory/Meteorological variables all sites.xlsx', sheet= 1, detectDates = TRUE)
df_met <- df_met[-c(1),]
df_met <- data.frame(Year = as.numeric(df_met$Year), DoY = as.numeric(df_met$DoY), 
                     Hour = as.numeric(df_met$Hour), LE = as.numeric(df_met$LE), 
                     Ta = as.numeric(df_met$Tair), VPD = as.numeric(df_met$VPD))
df_met$DateTime <- make_datetime(year = df_met$Year, min = round(df_met$Hour*60), tz='America/Chicago') + days(df_met$DoY-1)
df_met <- df_met[,!(names(df_met) %in% c('Year', 'DoY', 'Hour'))]
df_met <- df_met %>% select('DateTime', everything())
df_met[df_met==-9999] <- NA

# GPP data
df_GPP <- read.xlsx('C://Users/ppushpendra.STUDENT/Box/FVS_theory/30 min GPP wheat sites.xlsx', sheet=1)
df_GPP <- data.frame(Year = as.numeric(df_GPP$Year), DoY = as.numeric(df_GPP$DoY), 
                     Hour = as.numeric(df_GPP$Hour), GPP = as.numeric(df_GPP$GPP_C))
df_GPP$DateTime <- make_datetime(year = df_GPP$Year, min = round(df_GPP$Hour*60),  tz='America/Chicago') + days(df_GPP$DoY-1)
df_GPP <- df_GPP[,!(names(df_GPP) %in% c('Year', 'DoY', 'Hour'))]
df_GPP <- df_GPP %>% select('DateTime', everything())
df_GPP[df_GPP==-9999] <- NA

# Rn data
df_Rn <- read.xlsx('C://Users/ppushpendra.STUDENT/Box/FVS_theory/Rn data for wheat sites.xlsx', sheet=1, detectDates = TRUE)
df_Rn[df_Rn==-9999] <- NA
df_Rn$time <- chron::times(df_Rn$time)
df_Rn$DateTime <- as.POSIXct(paste(df_Rn$date, df_Rn$time), format="%Y-%m-%d %H:%M:%S",  tz='America/Chicago')
df_Rn <- data.frame(DateTime=df_Rn$DateTime, Rn=as.numeric(df_Rn$NR_Wm2_Avg))

# G fluxes Data
df_G <- read.xlsx('C://Users/ppushpendra.STUDENT/Box/FVS_theory/G fluxes wheat sites.xlsx', sheet=1, detectDates = TRUE)
df_G <- df_G[-c(1),]
df_G[df_G==-9999] <- NA
df_G$time <- chron::times(df_G$time)
df_G$DateTime <- as.POSIXct(paste(df_G$date, df_G$time), format="%Y-%m-%d %H:%M:%S",  tz='America/Chicago')
df_G <- data.frame(DateTime=df_G$DateTime, G=as.numeric(df_G$G))

# ELRE station data
df_ELRE <- read.csv('C://Users/ppushpendra.STUDENT/Box/FVS_theory/ELRE_Data_30min.csv')
df_ELRE[df_ELRE==-999] <- NA
df_ELRE$DateTime <- as.POSIXct(df_ELRE$TIME, format="%Y-%m-%dT%H:%M",  tz='America/Chicago')
df_ELRE <- data.frame(DateTime=df_ELRE$DateTime, Rain=(df_ELRE$RAIN*25.4), Pa = df_ELRE$PRES*3386.39/1000, SRad = df_ELRE$SRAD)
#
df <- left_join(df_met, df_GPP, by='DateTime') %>%
  left_join(., df_Rn, by='DateTime') %>%
  left_join(., df_G, by='DateTime') %>%
  left_join(., df_ELRE, by='DateTime')

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Add PET (using P-T 1972)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
df$LE_pot <- calc_PET_PT(Ta = df$Ta, Pa = df$Pa, Rn = df$Rn, G = df$G, alpha=1.26)$LE_pot
df$PET <- calc_PET_PT(Ta = df$Ta, Pa = df$Pa, Rn = df$Rn, G = df$G, alpha=1.26)$PET
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Data screening and quality control 
# --> (1) exclude defective entries and only daylight data (7 A.M. to 7 P.M.) 
# with positive Rn, GPP, ET, and VPD
# --> (2) exclude data from rainy days and several days that followed rainy days
# as follows: Two dry days following a rainy day will be excluded when P >= 2*PET, 
# otherwise one dry day will be excluded when P >= PET. Only the rainy days will be
# excluded when P < PET for the day.
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
df$Mask <- TRUE
df$Mask[(hour(df$DateTime) + minute(df$DateTime)/60) < 7 | (hour(df$DateTime) + minute(df$DateTime)/60) > 19 ] <- FALSE
df$Mask[df$LE < 0 | df$Rn < 0 | df$GPP < 0 | df$VPD < 0] <- FALSE

df$PET[!df$Mask] <- NA
df$Rain[!df$Mask] <- NA
df$PET[is.na(df$Rain)] <- NA
df$Rain[is.na(df$PET)] <- NA

df <- df %>% 
  group_by(Date = date(DateTime)) %>% 
  mutate(Rain_daily = mean(Rain, na.rm=T)) %>%
  mutate(PET_daily = mean(PET, na.rm=T))

df$Rain_daily <- df$Rain_daily*24*60/30    # 30 min-1 to day-1
df$PET_daily <- df$PET_daily*24*3600       # s-1 to day-1

df$Mask[df$Rain_daily > 0] <- FALSE

df <- df %>% as.data.table() # to simplify operations on df
df <- df %>% 
  mutate(Mask = case_when(Date %in% c(df[Rain_daily > PET_daily, Date],df[Rain_daily > PET_daily, Date]+1) ~ FALSE,
                          TRUE ~Mask ))
df <- df %>% 
  mutate(Mask = case_when(Date %in% c(df[Rain_daily > 2*PET_daily, Date],df[Rain_daily > 2*PET_daily, Date]+1, 
                                      df[Rain_daily > 2*PET_daily, Date]+2) ~ FALSE, TRUE ~Mask ))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
df <- subset(df, date(DateTime) >= "2016-11-01" & date(DateTime) <= "2017-05-31")
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
df$ET  <- LE2ET(df$LE, df$Ta)
df$GPP_mul_sqrt_VPD <- df$GPP*sqrt(df$VPD)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
df <- df[!(is.na(DateTime)),]
DateTime <- df$DateTime
Mask <- df$Mask
df[(df$Mask==FALSE), ] <- NA   # Select only screened data
df$DateTime <- DateTime
df$Date <- date(df$DateTime)
df$Mask <- Mask
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Visualization
plot(df$ET*24*3600, df$GPP_mul_sqrt_VPD*48, xlab = expression(ET ~ (kg ~ H[2]*O~m^-2 ~ day^-1)), ylab="", 
     cex.lab=1.5, cex.axis=1.2, xlim = c(0,20), ylim=c(0,250), pch=20)
title(ylab=expression(GPP%.%VPD^0.5~(gC~hPa^0.5~m^-2~day^-1)), line=2.2, cex.lab=1.5)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Fit 95th Quantile regression
#library(quantreg)
df$ET_30min <- df$ET*30*60   # mm s-1 to mm 30min-1
df$GPP_mul_sqrt_VPD_day <- df$GPP_mul_sqrt_VPD*48   # x 30min-1 to x day-1
#qt_fit <- rq(GPP_mul_sqrt_VPD_day ~ ET_day, tau = 0.95, data = df)
temp_df <- data.frame(x=df$ET_30min, y=df$GPP_mul_sqrt_VPD_day)
temp_df <- temp_df[complete.cases(temp_df), ]
qt_fit <- quantile_reg(temp_df$x, temp_df$y,  PolyDeg=0, tau=0.95, weights=NA)

#uWUEp <- qt_fit[[1]][2]
uWUEp <- qt_fit[[1]]
print(uWUEp)
x <- seq(0,16, 0.01)
y <- uWUEp*x #+ qt_fit[[1]][1]
par(new=TRUE)
lines(x,y,lty=2,col="red",lwd=3)
text(5, 240, expression(uWUE[p]== ~19.3125~gC%.%hPa^0.5/kg~H[2]*O), cex = 1.5)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# T:ET at half-hourly scale
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
temp <- df
temp$uWUEp <- uWUEp
temp$uWUEa <- temp$GPP_mul_sqrt_VPD_day/temp$ET_30min
temp$T_ET <- temp$uWUEa/temp$uWUEp
temp[temp==Inf] <- NA
temp$T_ET[temp$T_ET<=0] <- NA
temp$T_ET[temp$T_ET>=1] <- NA

# Let's compare with FVS outputs
df_FVS <- read.csv('C://Users/ppushpendra.STUDENT/Box/FVS_theory/DataPreprocessing_Raghav/Output/Site_1/const_ratio_2016_2018_all_vars.csv', skip = 1)
df_FVS <- data.frame(DateTime=df_FVS$X, LEt_FVS=df_FVS$LEt, LE_FVS=df_FVS$LE)
df_FVS$DateTime <- as.POSIXct(df_FVS$DateTime, format="%Y-%m-%d %H:%M:%S", tz='America/Chicago')
df_FVS$T_ET_FVS <- df_FVS$LEt_FVS/df_FVS$LE_FVS

temp <- data.frame(DateTime = temp$DateTime, T_ET_uWUE=temp$T_ET)
temp <- left_join(df_FVS, temp, by='DateTime')
plot(temp$T_ET_FVS, temp$T_ET_uWUE, xlab = "T:ET (FVS)", ylab="", 
     cex.lab=1.5, cex.axis=1.2, xlim = c(0,1), ylim=c(0,1), pch=20)
title(ylab="T:ET (uWUE)", line=2.2, cex.lab=1.5)
title(main="half-hourly estimates")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# T:ET at daily scale
# uWUEa will be estimated as the linear regression slope from a moving window spanning one day
# T:ET will be estimated only for days when there are at least 10 effective entries
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
temp <- df
temp <- temp %>% as.data.table() # to simplify operations on df
temp <- temp %>% 
  group_by(Date) %>% 
  mutate(ET_30min_vals_count = sum(!(is.na(ET_30min)))) %>%  # Count number of effective entries of ET data
  mutate(GPP_mul_sqrt_VPD_day_count = sum(!(is.na(GPP_mul_sqrt_VPD_day)))) # Count number of effective entries of GPP*sqrt(VPD) data

temp$Mask[(temp$ET_30min_vals_count < 10 | temp$GPP_mul_sqrt_VPD_day_count < 10)] <- FALSE
DateTime <- temp$DateTime
Mask <- temp$Mask
temp[(temp$Mask==FALSE), ] <- NA   # Select only screened data
temp$DateTime <- DateTime
temp$Date <- date(temp$DateTime)
temp$Mask <- Mask
#-------------------------------------------------------------------------------
# Linear Regression for uWUEa for each day
temp$uWUEa <- NA
for(i in 1:length(unique(temp$Date))){
  ET_30min <- temp$ET_30min[temp$Date==unique(temp$Date)[i]]
  GPP_mul_sqrt_VPD_day <- temp$GPP_mul_sqrt_VPD_day[temp$Date==unique(temp$Date)[i]]
  if(sum(!(is.na(ET_30min)))>1){
    temp.df <- data.frame(x=ET_30min, y=GPP_mul_sqrt_VPD_day)
    temp.df <- temp.df[complete.cases(temp.df),]
    linearMod <- lm(y ~ x, data=temp.df)
    uWUEa <- linearMod["coefficients"][[1]][2]
    temp$uWUEa[temp$Date==unique(temp$Date)[i]] <- uWUEa
  }
}
#-------------------------------------------------------------------------------
temp$uWUEp <- uWUEp
temp$T_ET <- temp$uWUEa/temp$uWUEp
# Let's compare with FVS outputs
df_FVS <- read.csv('C://Users/ppushpendra.STUDENT/Box/FVS_theory/DataPreprocessing_Raghav/Output/Site_1/const_ratio_2016_2018_all_vars.csv', skip = 1)
df_FVS <- data.frame(DateTime=df_FVS$X, LEt_FVS=df_FVS$LEt, LE_FVS=df_FVS$LE)
df_FVS$DateTime <- as.POSIXct(df_FVS$DateTime, format="%Y-%m-%d %H:%M:%S", tz='America/Chicago')
df_FVS$LEt_FVS[((hour(df$DateTime) + minute(df$DateTime)/60) < 7 | (hour(df$DateTime) + minute(df$DateTime)/60) > 19)] <- NA
df_FVS$LE_FVS[((hour(df$DateTime) + minute(df$DateTime)/60) < 7 | (hour(df$DateTime) + minute(df$DateTime)/60) > 19)] <- NA
df_FVS <- df_FVS %>% 
  group_by(Date=date(DateTime)) %>% 
  mutate(vals_count = sum(!(is.na(LEt_FVS)))) # Count number of effective entries of FVS output between 7AM-7PM

df_FVS$LEt_FVS[(df_FVS$vals_count < 10)] <- NA
df_FVS$LE_FVS[(df_FVS$vals_count < 10)] <- NA
df_FVS <- df_FVS %>% group_by(Date) %>%
  summarise_each(funs(mean(., na.rm=T)))
df_FVS <- data.frame(Date=df_FVS$Date, T_ET_FVS=df_FVS$LEt_FVS/df_FVS$LE_FVS)

temp <- data.frame(DateTime = temp$DateTime, T_ET_uWUE=temp$T_ET)
temp <- temp %>% group_by(Date=date(DateTime)) %>%
  summarise_each(funs(mean(., na.rm=T)))
temp <- data.frame(Date=temp$Date, T_ET_uWUE=temp$T_ET_uWUE)

temp <- left_join(temp, df_FVS, by='Date')
plot(temp$T_ET_FVS, temp$T_ET_uWUE)



df <- read.csv('C://Users/ppushpendra.STUDENT/Box/test.csv')
x <- df$x
y <- df$y
data <- data.frame(x=x, y=y)
linearMod <- lm(y~x, data=data)
uWUEa <- linearMod["coefficients"][[1]][2]
