
############## Setup, Graph plotting and basic identification #############

# remove existing environment variables
rm(list=ls())

# set working directory
setwd("/Users/davidhui/Documents/GitHub/uon/MATH4022 Time Series and Forecasting/")

# import the a23_nox.csv
raw_data = read.csv("materials/a23_nox.csv")
head(raw_data)  # see the first few rows
n = nrow(raw_data)  # number of observations
n

# convert into ts object
data = ts(raw_data$daily_mean_nox)

mean(data)  # mean of the data is 65.26074

# plot the row ts plot
ts.plot(data, main="Daily mean nitrous oxides level from 1/2/2017 to 30/9/2017",
        xlab="Time",
        ylab=expression(paste("mean nitrous oxides level (in ug/",m^3,")")))

# from the ts plot, there does not appear to have trends and seasonal components; 
# also, the plotted data appears to have constant mean and constant variance

par(mfrow=c(2,1))
acf(data, main="sample ACF plot", ylim=c(-1,1), lag.max=30)
pacf(data, main="sample PACF plot", ylim=c(-1,1), lag.max=30)
par(mfrow=c(1,1))

# sample ACF and sample PACF decay to zero gradually as lags increase; 
# it seems plausible to claim that the time series is stationary.

# to identify the order of the ARMA(p,q) process

# from the sample ACF plot, there do not have obvious cut-off after some lag
# > so the data appears to include AR components
# from the sample PACF plot, there "seems" a drop in PACF value after lag 1, 
# > however, the sample PACF at lag 2 is still slightly lie outside +- 2 /sqrt(n) bound,
# > so this is in doubt to claim the PACF show obvious cut-off.

# Initially, I will suggest and AR(1) first
# and then try to introduce different ARMA model


# so let {X_t} be the original process
# note that the mean of the process is non-zero, 
# so in our arima(...) model fitting, the estimate of the intercept will be non-zero

########### Identify the order of an ARMA(p,q) ###########

#### first take a look at my initial suggested model AR(1)
# (X_t - mu) = phi1 * (X_t-1 - mu) + Z_t ; var(Z_t) = sigma_z^2
arima(data, order=c(1,0,0), method="ML")
# by the hypothesis: H0: phi1 = 0 vs H1: phi1 != 0
# since the estimate = 0.4994 with s.e. = 0.0555 
#   => Test Statistic = abs(0.4994/0.0555) = 8.9981 > 2
# we can see that the test statistic for phi1 (ar1) are greater than 2
# and we should keep phi1 into the model

## since replicate these steps are quite time consuming; 
## so I am going to write a function to print the fit summary 
##    for different AR and MA order
## data = the given time series object
## p = order of AR components
## q = order of MA components
ARMA_fit = function(data, p, q){
  fit = arima(data, order=c(p,0,q), method="ML")
  print("===== fit summary =====")
  print(fit)
  print("===== test statistic =====")
  # coef(fit) returns the parameter estimation; 
  # sqrt(diag(vcov(fit))) returns the s.e. of each parameter est.
  test_statistic = (coef(fit) / sqrt(diag(vcov(fit))))
  # no need to print "intercept
  print(test_statistic[1:(p+q)])
  
  # return the arima fit object, later I will need to extract the AIC
  return (fit)
}

# to see if this works, repeat for AR(1)
# this print the same results
fit_AR1 = ARMA_fit(data, 1, 0)

## try fitting different orders, each time introduce one more parameter
## because of the principle of parsimony
## we shouldn't add too much parameters to the model
## since simple hypothesis is to test the significance of the newly introduced parameter, 
##  > given the others are presented
## so I keep trying until I found the newly added parameter is non-significant
## NOTE that this may comes from the problem of multiple testing; 
## nevertheless, this should be helpful for me to find "possible" model first
## later I will use AIC to check if 
## introducing more parameters will have significant improvement

# AR(2)
fit_AR2 = ARMA_fit(data, 2, 0)  
# phi2 is significant, given phi1 is presented

## AR(3)
fit_AR3 = ARMA_fit(data, 3, 0)  
# phi3 is not significant given phi1 and phi2 are presented

# so for pure AR process, I will end up with AR(2)


# next start to fit ARMA
fit_ARMA11 = ARMA_fit(data, 1, 1) 
# theta1 is significant, given phi1 is presented

fit_ARMA21 = ARMA_fit(data, 2, 1) 
# phi2 is not significant given phi1 and theta1 are presented

fit_ARMA12 = ARMA_fit(data, 1, 2)
# theta2 is not significant given phi1 and theta1 are presented

fit_ARMA22 = ARMA_fit(data, 2, 2) 
# throws error
# should because of the parameter estimate fail 
# to fall into "appropriate region" for stationarity/invertibility
# so for ARMA process, I will end up with ARMA(1,1)

# try to fit just MA as well
fit_MA1 = ARMA_fit(data, 0, 1)
# theta1 is significant

fit_MA2 = ARMA_fit(data, 0, 2) 
# theta2 is significant, given theta1 is presented

fit_MA3 = ARMA_fit(data, 0, 3) 
# theta3 is not significant given theta1 and theta2 are presented

# so for pure MA process, I will end up with MA(2)

## combining the above information, we can just compare the AIC 
## of AR(1); AR(2); MA(1); MA(2); ARMA(1,1)
## since the sample size is 242; which is not very small; 
## so I will just use AIC instead of corrected AIC
aic_df = data.frame(
  model=c("AR(1)", "AR(2)", "MA(1)", "MA(2)", "ARMA(1,1)"),
  AIC = c(AIC(fit_AR1), AIC(fit_AR2), AIC(fit_MA1), AIC(fit_MA2), AIC(fit_ARMA11))
)

aic_df

# since for model selection procedure; the small AIC indicates a "better" model
# I start with AR(1), where AIC = 2200.872
# MA(1) gives a larger AIC value, so I would not take MA(1) into account, 
# and hence MA(2) is not reasonable to consider
# to compare AR(1) with AR(2) and ARMA(1,1), AR(2) gives the smallest AIC. 
# although the difference is not very large (2200.872 for AR(1) vs 2197.361 for AR(2))
# however, since phi2 is significant at 5% significant level
# so I think AR(2) is a reasonable model for this dataset

# extract all parameter estimate from fit_AR2
phi1 = signif(coef(fit_AR2)[[1]],4)
phi2 = signif(coef(fit_AR2)[[2]],4)
mu = signif(coef(fit_AR2)[[3]],4)
sigma_z_square = signif(fit_AR2$sigma2,4)

# so the suggested model is 
# (X_t-65.33)=0.5753*(X_t-1-65.33)-0.1501*(X_t-2-65.33)+Z_t;
# with var(Z_t) = 496.6

############## Model Diagnostic #############
# ASSUMING OUR MODEL IS CORRECT

# finally for diagnostic, plot the following:
# ACF plot: sample ACF with overlay theoretical ACF 
# PACF plot: sample PACF with overlay theoretical PACF
# residual ts plot
# residual acf plot
# Ljung-Box test with plot
source("materials/Computer Practical 3 - R function for Ljung-Box test.R")
par(mfrow=c(2,1))

# sample ACF plot
acf(data, main="sample ACF plot", ylim=c(-1,1), lag.max=30)
# overlay theoretical ACF, using the parameter estimate
lines(0:30, ARMAacf(ar = c(phi1, phi2), ma = c(0), lag.max = 30), col="red")
# legend on top right
legend(x="topright", 
       legend=c("sample ACF", "theoretical ACF"),
       col=c("black", "red"),
       lty=1)

# sample PACF plot
pacf(data, main="sample PACF plot", ylim=c(-1,1), lag.max=30)
# overlay theoretical PACF, using the parameter estimate
lines(1:30, ARMAacf(ar = c(phi1, phi2), ma = c(0), lag.max=30, pacf=TRUE), col="red")
# legend on top right
legend(x="topright", 
       legend=c("sample PACF", "theoretical PACF"),
       col=c("black", "red"),
       lty=1)

par(mfrow=c(1,1))

# ACF and PACF are fitting quite good using our fitted model



## next see the residual plot and Ljung-Box test plot
resid_fit_AR2 = residuals(fit_AR2)

## using the material from Computer practical 3 for Ljung-Box test
# order of AR: p = 2
# order of MA: q = 0
# max.k = 12 so test degrees of freedom up to 10 (10 + p + q = 12)
LB_fit_AR2 = LB_test(resid_fit_AR2, max.k = 12, p = 2, q = 0)

par(mfrow=c(3,1))

# residuals time plot
plot(resid_fit_AR2, main="Residual Plot for the AR(2) fit", xlab="Time", ylab="Residual")

# residuals ACF plot
acf(resid_fit_AR2, main="ACF plot of residual for AR(2) fit", ylim=c(-1,1))

# LB-test result plotting
plot(LB_p_value ~ deg_freedom, data = LB_fit_AR2, 
     main="Ljung-Box test: p-value vs degrees of freedom",
     xlab="degrees of freedom", ylab="p-value")

# add a reference line (horizontal 0.05)
abline(h=0.05, col="blue", lty=2)
par(mfrow=c(1,1))



# from the sample ACF plot, at lag 7, 14, 21
# there are ACF value lie outside the bound, but in fact they are close to the bound
# so it is insufficient to claim the residuals appears to have "seasonal" components
# from the Ljung-Box test plot, 
# Because the LB test is used to test 
# H0: residuals are independent vs H1: residuals are dependent
# we expect the p-value to be above 0.05 to retain H0
# However, from this model fit, p-values are below 0.05 on and after df=5
# it indicates the residuals are dependent


## next I define a helper function to plot these diagnostic plots
## then see the others fitted model

## data: the time series object
## fitted_model: the arima model
## p: order of AR component
## q: order of MA component
resid_diagnostic_plot = function(data, fitted_model, p, q){
  # extract the residuals of the model fit
  resid = residuals(fitted_model)
  
  # get the Ljung-Box test p-value for each degrees of freedom
  LB_fit = LB_test(resid, max.k = 10 + p + q, p = p, q = q)

  par(mfrow=c(3,1))
  
  # plot residuals plot
  plot(resid, main=paste("Residual Plot for the ARMA(",p,",",q,") fit"), 
       xlab="Time", ylab="Residual")

  # plot residuals' ACF plot
  acf(resid, main=paste("ACF plot of residual for ARMA(",p,",",q,") fit"), 
      ylim=c(-1,1))
  
  # plot Ljung-Box test results: p-value against df
  plot(LB_p_value ~ deg_freedom, data = LB_fit, 
       main="Ljung-Box test: p-value vs degrees of freedom",
       xlab="degrees of freedom", ylab="p-value") 
  # add a reference line: at 0.05 significance level
  abline(h=0.05, col="blue", lty=2)
  
  par(mfrow=c(1,1))
  
}

## my proposed model AR(2)
resid_diagnostic_plot(data, fit_AR2, p=2, q=0)  

## rest of the fitted model
resid_diagnostic_plot(data, fit_AR1, p=1, q=0)
resid_diagnostic_plot(data, fit_AR3, p=3, q=0)
resid_diagnostic_plot(data, fit_MA1, p=0, q=1)
resid_diagnostic_plot(data, fit_MA2, p=0, q=2)
resid_diagnostic_plot(data, fit_ARMA11, p=1, q=1)
resid_diagnostic_plot(data, fit_ARMA12, p=1, q=2)
resid_diagnostic_plot(data, fit_ARMA21, p=2, q=1)

### none of these satisfied Ljung-Box test criteria
### so I may just end up with propose AR(2) to be the best model here

