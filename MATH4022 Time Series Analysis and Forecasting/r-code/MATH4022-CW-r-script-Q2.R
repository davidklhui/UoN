
######## Q2

###################### BASIC INVESTIGATION ########

rm(list=ls()) # remove environment variables

# set working directory
setwd("/Users/davidhui/Documents/GitHub/uon/MATH4022 Time Series and Forecasting/")

# read the eng_car_reg.csv
raw_data = read.csv("materials/eng_car_reg.csv")
head(raw_data)  # see the first few rows
n = nrow(raw_data)  # number of observations
n # n = 87

# convert the imported data into time series object 
# frequency = 4 (quartly)
# with start at 1st quarter of 2001

data = ts(raw_data$no_new_regs, start=c(2001,1), frequency = 4)
data

# time series plot of the data
ts.plot(data, main="Quarterly number of new car registrations in England (in Thousands)",
        xlab = "Year",
        ylab = "numbers of new car registrations (in thousands)")

# it shows a seasonal pattern: 
#   initially Q1 would be high, 
#   then in Q2 would be a slight drop,
#   then in Q3 would rise again,
#   finally in Q4 would be a larger drop and reach the minimum of the year
# and repeats this pattern year by year generally

# it appears to have obvious outliers (in 2008 Q3, and 2020 Q2)
# in these 2 quarter, it has a large drop
# possible reasons: because 2008 Q3 had financial crisis; 
# and 2020 Q2 had COVID-19 pandemic
# these events may have a large impact on new car registrations.

# except for where the outliers occurred,
# the variability are generally the same 
# (no systematic change in variability)
# the mean is changing over different period. 
# from 2001 to 2010, it appears to have a decreasing trend linearly
# from 2010 to 2017, it appears to have a increasing trend linearly
# finally from 2017 and afterwards, it showed a decreasing trend again

# to see the ACF and PACF plot
par(mfrow=c(2,1))
acf(data, main="ACF plot", lag.max=100)
pacf(data, main="PACF plot", lag.max=100)
par(mfrow=c(1,1))
# from the ACF plot, we can see that there are spikes for each lag multiples of 2
# to best describe this phenomenon, 
# this is because from the data shows seasonal pattern
# so for every even quarters would have a drop; 
# and for every odd quarters would have a rise
# so this will give a strong correlation for each lag multiples of 2
# however, we can see that the lag multiples of 4 would be a bit "larger"
# because of the influence of the seasonal period = 4


#################### MODEL IDENTIFICATION #####################


#### to start guessing the order of the SARIMA model ARIMA(p,0,q)x(P,1,Q)[4]
#### to start, because there are obvious seasonal trends, first remove seasonal pattern
data_diff4 = diff(data, lag=4)
layout(matrix(c(1,1,2,3), nc=2, byrow=T))
ts.plot(data_diff4, main="Time Series plot of seasonal differenced data (period=4)")
acf(data_diff4, ylim=c(-1,1), main="ACF plot of seasonal differenced data", lag.max=30)
pacf(data_diff4, ylim=c(-1,1), main="PACF plot of seasonal differenced data", lag.max=30)
par(mfrow=c(1,1))

# first from the time-series plot, even though there are some 
# rapid fluctuation around 2010 and 2020,
# however apart from those that I determined as "outliers", 
# we see that the rest parts of the plot looks "stationary"
# from the sample ACF and PACF plot, 
# because both declines to 0 as lag increases

# we can see that both ACF and PACF plot, have a clear spike at lag 4 (seasonal period)
# except for the clear spike at lag 4, 
# both ACF and PACF plot doesn't seem to have obvious cut-off
# it might show evidence that both AR and MA component exists in the model



################## MODEL FITTING ##################

## data = raw data
## p = non-seasonal AR order
## d = non-seasonal difference order
## q = non-seasonal MA order
## P = seasonal AR order
## D = seasonal difference order
## Q = seasonal MA order
## S = period
SARIMA_fit = function(data, p, d, q, P, D, Q, S){
  fit = arima(data, order=c(p,d,q), seasonal = list(order=c(P,D,Q), period=S),
              method="ML")
  print("===== fit summary =====")
  print(fit)
  print("===== test statistic =====")
  # coef(fit) returns the parameter estimation; 
  # sqrt(diag(vcov(fit))) returns the s.e. of each parameter est.
  test_statistic = (coef(fit) / sqrt(diag(vcov(fit))))
  # no need to print "intercept
  print(test_statistic)
  
  # return the arima fit object, later I will need to extract the AIC
  return (fit)
}


df = as.numeric()

## I only try up to order 2
## and because I believe that the period is S = 4, 
## only need perform 1st order of seasonal diff D = 1
## no need to further perform non-seasonal diff d = 0
## because arima(...) may though error if parameter out of specific range
## so there is another tryCatch block to handle this auto-search method
for(p in 0:2){
  for(q in 0:2){
    for(P in 0:2){
      for(Q in 0:2){
        print(paste("fitting: p=",p,"q=",q,",P=", P, ",Q=", Q))
        tryCatch({
          fit = SARIMA_fit(data, p, 0, q, P, 1, Q, 4)
          
          # check if this auto-search provides any error or warning or NaN
          # if so , ignore them
          test_statistic = (coef(fit) / sqrt(diag(vcov(fit))))
          if(sum(is.na(test_statistic)) == 0){
            
            ## we need the nonseasonal orders (p,q), seasonal orders (P,Q), 
            tmp = data.frame(p = p, q = q, P = P, Q= Q, AIC = AIC(fit))
            df = rbind(df, tmp)
          }
        }, error = function(err){
          message(err)
          print(paste("error in fitting: p=",p,"q=",q,",P=", P, ",Q=", Q))
        }, warning = function(warn){
          message(warn)
          print(paste("warning in fitting: p=",p,"q=",q,",P=", P, ",Q=", Q))
        })
      }
    }
  }
}

# see the result
df

library("dplyr")

# sort by AIC
df %>% arrange(AIC)
# from this greedy auto-search, we can find that 
# there are 76 out of 81 iterations give valid AIC
# we can checkout the rest of them

# see the first 10 only
df %>% arrange(AIC) %>% head(10)


# choose ARIMA(1,0,1)(0,1,1)[4]
fitted_model = SARIMA_fit(data, p=1, d=0, q=1, P=0, D=1, Q=1, S=4)


################## MODEL DIAGNOSTIC ####################
### 1. ACF and PACF with overlay theoretical ACF/PACF (on seasonal diff data)
### 2. residuals plot
### 3. residuals ACF plot
### 4. Ljung-Box test


## 1. ACF and PACF with overlay theoretical ACF/PACF
# using our model ARIMA(1,0,1)(0,1,1)[4]
# can be read as ARMA(1,5) for the seasonal differenced data

phi1 = 0.9845;
theta1 = -0.5401;
Theta1 = -0.9360;
S = 4;
ar = c(phi1)
ma = c(theta1, 0, 0, Theta1, theta1*Theta1)

# sample ACF plot on the seasonal differenced data
acf(data_diff4, main="sample ACF plot of seasonal diff data", 
    ylim=c(-1,1), lag.max=30)
# overlay theoretical ACF, using the parameter estimate
# note that we need the "Lag" divided by the seasonal period
lines((0:30)/S, ARMAacf(ar = ar, ma=ma, lag.max = 30), col="red")
# legend on top right
legend(x="topright", 
       legend=c("sample ACF", "theoretical ACF"),
       col=c("black", "red"),
       lty=1)

# sample PACF plot on the seasonal differenced data
pacf(data_diff4, main="sample PACF plot of seasonal diff data"
     , ylim=c(-1,1), lag.max=30)
# overlay theoretical PACF, using the parameter estimate
lines((1:30)/S, ARMAacf(ar = ar, ma = ma, lag.max = 30, pacf=TRUE)
      , col="red")
# legend on top right
legend(x="topright", 
       legend=c("sample PACF", "theoretical PACF"),
       col=c("black", "red"),
       lty=1)

# ACF and PACF are fitting quite good using our fitted model

### residuals plot
fitted_model_resid = residuals(fitted_model)
ts.plot(fitted_model_resid, 
        main="Residuals plot")

### residuals ACF plot
acf(fitted_model_resid, 
    main="Residuals ACF plot")


### Ljung-Box test
source("materials/Computer Practical 3 - R function for Ljung-Box test, SARIMA MODEL.R")
## using the material from Computer practical 3 for Ljung-Box test
# non-seasonal order of AR: p = 1
# non-seasonal order of MA: q = 1
# seasonal order of AR: P = 0
# seasonal order of MA: Q = 1
# max.k = 13 so test degrees of freedom up to 10 (10 + p + q + P + Q= 13)
LB_fitted_model_resid = LB_test_SARIMA(fitted_model_resid, max.k=13, p=1, q=1, P=0, Q=1)

# plot the Ljung-Box test resulting plot
plot(LB_p_value ~ deg_freedom, data = LB_fitted_model_resid, 
     main="Ljung-Box test: p-value vs degrees of freedom",
     xlab="degrees of freedom", ylab="p-value",
     ylim=c(0,1))

# add a reference line (horizontal 0.05)
abline(h=0.05, col="blue", lty=2)


################# FORECASTING #################

library("forecast")
## using forecast method
## data = raw data
## model = using our ARIMA(1,0,1)(0,1,1)[4]
## h = 4, our required future time step for forecasting
## level = 95, the confidence level
forecast_value = forecast(data, model=fitted_model, h=4, level=95)


# check out the forecasted value with uncertainty
# uncertainty = (upper95% bound - lower95% bound ) / 2
fc = forecast_value %>% 
  as.data.frame(.) %>%
  rename("Lower" = "Lo 95", "Higher" = "Hi 95") %>%
  mutate("Uncertainty" = (Higher - Lower)/2)

fc


# combine raw data with forecasting mean and uncertainty
data_with_fc = cbind(forecast_value$mean, 
                     forecast_value$upper, 
                     forecast_value$lower, 
                     data)

# plot those data
# together with associated uncertainty
autoplot(data_with_fc,
         main="Quarterly number of new cars registered with prediction",
         xlab="Year", ylab="number of new cars registeration (in thousands)")



