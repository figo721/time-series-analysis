###################################################################################
### Advanced Time Series Analysis ######
### First Name: Jaechan 
### Last  Name: Park
### please contact at figo721@hotmail.com regarding the source code 


rm(list = ls())

library("quantmod")
library(CADFtest)
library(forecast)
library(fGarch)
library(vars)
library(urca)


nasdaq<-getSymbols('^IXIC', src='yahoo', auto.assign = FALSE)
dfnasdaq<-as.data.frame(nasdaq$IXIC.Close)

amazon_data<-getSymbols('AMZN', src='yahoo',auto.assign = FALSE)
dfamazon<-as.data.frame(amazon_data$AMZN.Close)

nasdaq_ts<- ts(dfnasdaq, frequency = 365,start = c(2010,1))
summary(nasdaq_ts)
plot.ts(nasdaq_ts)  # Has upward Trend , no seasonality, random walk with drift,  


##################################################################################
##### Difference of nasdaq.close index ####
dnasdaq_ts<-diff(nasdaq_ts)
ts.plot(dnasdaq_ts)
acf(dnasdaq_ts)
pacf(dnasdaq_ts)

x<-dim(nasdaq)
TREND<-1:x  # has to be modified the end value of the "TREND" since the data set updates simultaneously
fit<-lm(dnasdaq_ts~TREND)
summary(fit)
ts.plot(fit$residuals)  # Heteroscedasticity

max.lag<-round(sqrt(length(dnasdaq_ts)))
mydata<-CADFtest(dnasdaq_ts, type="drift", criterion = "BIC",max.lag.y=max.lag)
summary(mydata)
# P-value of the test is 0.00<5% => we reject the null hypothesis and conclude that there is a stochastic trend.

# to check white noise for Ljung-Box test.
Box.test(dnasdaq_ts, lag=max.lag, type="Ljung-Box")
# since p-value=0.00<5% => reject the null hypothesis and conclude that it is not a white noise.


#################################################################################



###################################################################################
##### Log transformation ###
lognasdaq_ts<-log(nasdaq_ts)   
ts.plot(lognasdaq_ts)    
y<-dim(lognasdaq_ts)
TREND<-1:y # has to be modified the end value of the "TREND" since the data set updates simultaneously
fit1<-lm(lognasdaq_ts~TREND)
summary(fit1)
ts.plot(fit1$residuals)  # Heteroscedasticity

acf(lognasdaq_ts)
pacf(lognasdaq_ts)

##Unit root test to check for its stationary
max.lag<-round(sqrt(length(lognasdaq_ts)))
CADFtest(lognasdaq_ts, type="drift", criterion = "BIC",max.lag.y=max.lag)
# p-value =0.9214 > 0.05 -> we do not reject the null hypothesis and conclude that it is a stochastic trend.

# to check white noise for Ljung-Box test.
Box.test(lognasdaq_ts, lag=max.lag, type="Ljung-Box")
# since p-value=0.00<5%, it is not a white noise.

##################################################################################



#################################################################################
### log difference
dlognasdaq_ts<-diff(log(nasdaq_ts))
ts.plot(dlognasdaq_ts)

# Checking white noise with acf and partial acf plots. 
acf(dlognasdaq_ts)
pacf(dlognasdaq_ts, lag.max = max.lag)
# this is not a white noise plot result.

z<-dim(dlognasdaq_ts)
TREND<-1:z  # has to be modified the end value of the "TREND" since the data set updates simultaneously
fit2<-lm(dlognasdaq_ts~TREND)
summary(fit2)
ts.plot(fit2$residuals)  # seems to be constant variance

#Unit root test to check for its stationary
max.lag<-round(sqrt(length(nasdaq_ts)))
CADFtest(dlognasdaq_ts, type="drift", criterion = "BIC",max.lag.y=max.lag)
# P-value of the test is 0.00<5% => we can conclude that it is stationary.

# to check white noise for Ljung-Box test.
Box.test(dlognasdaq_ts, lag=max.lag, type="Ljung-Box")
# since p-value=0.00<5%, it is not a white noise.

###################################################################################



###################################################################################
##############################
### ARIMA(1,1,0)
fit_arma1<-arima(lognasdaq_ts, order=c(1,1,0))
fit_arma1

max.lag1<-round(sqrt(length(fit_arma1)))
plot(fit_arma1$residuals)
acf(fit_arma1$residuals)
pacf(fit_arma1$residuals)
Box.test(fit_arma1$residuals, lag=max.lag, type="Ljung-Box")
qqnorm(fit_arma1$residuals)
qqline(fit_arma1$residuals)



#############################################################################
## ARIMA (2,1,0) model

fit_arma2<-arima(lognasdaq_ts, order=c(2,1,0))
fit_arma2
summary(fit_arma2)

max.lag<-round(sqrt(length(fit_arma2)))
plot(fit_arma2$residuals)
acf(fit_arma2$residuals)
pacf(fit_arma2$residuals)
Box.test(fit_arma2$residuals, lag=max.lag, type="Ljung-Box")
qqnorm(fit_arma2$residuals)
qqline(fit_arma2$residuals)

## Compare the model, ARIMA(2,1,0) and ARIMA (1,1,0)

AIC(fit_arma1)
AIC(fit_arma2)
AIC(fit_arma1, k=max.lag)
AIC(fit_arma2, k=max.lag1)




#################################################################
### Forecast with ARIMA (2,1,0)  ####

myforecast_arma<-predict(fit_arma2, n.ahead = 365)
expected<-myforecast_arma$pred

lower<-myforecast_arma$pred-qnorm(0.975)*myforecast_arma$se
upper<-myforecast_arma$pred+qnorm(0.975)*myforecast_arma$se
cbind(lower,expected,upper)

plot.ts(lognasdaq_ts,xlim=c(2010,2020))
lines(expected,col="red")
lines(lower, col="blue")
lines(upper, col="blue")






# Compare the forecasts for the ARIMA (2,1,0) model and the ARIMA (1,1,0 model)
# by usin Diebold-mariano test statistics 
# the test has to be using the newey-west standard error which is wider.

y<-lognasdaq_ts
S=round(0.75*length(y))
h=1  # one step ahead
error1.h<-c()
for (i in S:(length(y)-h))
{
  mymodel.sub<-arima(y[1:i], order = c(2,1,0),seasonal=c(0,0,0))
  predict.h<-predict(mymodel.sub,n.ahead=h)$pred[h]
  error1.h<-c(error1.h,y[i+h]-predict.h)
}

error2.h<-c()
for (i in S:(length(y)-h))
{
  mymodel.sub<-arima(y[1:i], order = c(1,1,0),seasonal=c(0,0,0))
  predict.h<-predict(mymodel.sub,n.ahead=h)$pred[h]
  error2.h<-c(error2.h,y[i+h]-predict.h)
}
cbind(error1.h,error2.h)


MAE1<-mean(abs(error1.h))
MAE2<-mean(abs(error2.h))


dm.test(error1.h,error2.h,h=h,power=1)
## Since the p-value =0.02816, we can conclude that the two ARIMA model seems to be different.
## ARIMA (2,1,0) will be use for the forecasting.

##################################################################################################

## GARCH model
### arima(2,1,0)-garch(1,1) ###


fit_garch1<-garchFit(~arma(2,0)+garch(1,1),data=dlognasdaq_ts)
summary(fit_garch1)
plot(fit_garch1)


fit_garch2<-garchFit(~arma(2,0)+garch(1,1),data=dlognasdaq_ts, cond.dist = "QMLE")
summary(fit_garch2)
plot(fit_garch2)



myforecast_arma1<-predict(fit_garch2, n.ahead = 365)
expected<-myforecast_arma1$pred

lower<-myforecast_arma1$pred-qnorm(0.975)*myforecast_arma1$se
upper<-myforecast_arma1$pred+qnorm(0.975)*myforecast_arma1$se
cbind(lower,expected,upper)

plot.ts(lognasdaq_ts,xlim=c(2010,2020))
lines(expected,col="red")
lines(lower, col="blue")
lines(upper, col="blue")


########################################################################################

## Multivariate  ####
## Dynamic model#####

amazon_ts<- ts(dfamazon, frequency = 365,start = c(2010,1))
summary(amazon_ts)
plot.ts(amazon_ts)  # Has upward Trend , no seasonality, random walk with drift,  

logamazon_ts<-log(amazon_ts)
plot.ts(logamazon_ts) # non-stationary


# Show multiple time series graph(plots)
ts.plot(lognasdaq_ts, logamazon_ts, col=c("black","blue"))



dlogamazon_ts<-diff(log(amazon_ts))
plot.ts(dlogamazon_ts)

# multiple time series plot for difference 
ts.plot(dlognasdaq_ts, dlogamazon_ts, col=c("black", "blue"))



# Unit-root test 
max.lag<-round(sqrt(logamazon_ts))
CADFtest(logamazon_ts,type="drift",criterion="BIC",max.lag.y=max.lag)
# Do not reject the null hypothesis and it is non-stationary time series.


CADFtest(dlogamazon_ts,type="drift",criterion="BIC",max.lag.y=max.lag)
# reject the null hypothesis and concoldue that it is a stationary time series.
# first order difference of log-transformed time series is a stationary. 


# Estimate a linear regression
fit_diml<-lm(lognasdaq_ts~logamazon_ts)
summary(fit_diml)
plot.ts(fit_diml$residuals)  # looks like a stationary.
acf(fit_diml$residuals)
pacf(fit_diml$residuals)
max.lag<-round(sqrt(length(nasdaq_ts)))
Box.test(fit_diml$residuals,lag=max.lag, type = "Ljung-Box")
# since the result of the Q-test is that p-value = 0.00 which means rejecting the null
# hypothesis and the model is not valid.


# Then we will estimate a Autoregressive Distributed Lag model of order 3 (ADLM(3)) and check its validity.

lag <- 3
n <- length(lognasdaq_ts)
lognasdaq.0 <- lognasdaq_ts[(lag+1):n]
logamazon.0<-logamazon_ts[(lag+1):n]
lognasdaq.1 <- lognasdaq_ts[lag:(n-1)]
logamazon.1<-logamazon_ts[lag:(n-1)]
lognasdaq.2 <- lognasdaq_ts[(lag-1):(n-2)]
logamazon.2<-logamazon_ts[(lag-1):(n-2)]
lognasdaq.3 <- lognasdaq_ts[(lag-2):(n-3)]
logamazon.3<-logamazon_ts[(lag-2):(n-3)]

fit_adlm<-lm(lognasdaq.0~lognasdaq.1+lognasdaq.2+lognasdaq.3+logamazon.0+logamazon.1+logamazon.2+logamazon.3)
plot.ts(fit_adlm$residuals)
acf(fit_adlm$residuals)
Box.test(fit_adlm$residuals,lag=max.lag, type = "Ljung-Box")
# since p-value=0.000, we reject the null hypothesis which is assuming a white noise.


# Test for Granger causality

fit_adlm_nox<-lm(lognasdaq.0~lognasdaq.1+lognasdaq.2+lognasdaq.3)
gran_cal<-anova(fit_adlm, fit_adlm_nox)
gran_cal

# We obtain a p-value = 0.00 < 0.05
# thus we do reject H0 of no Granger Causality. 
# We conclude that logamazon has incremental explanatory power in predicting dlognasdaq.




## Testing whether lognasdaq and logamazon are cointegrated using the Engle-Granger test
fit_reg<-lm(lognasdaq_ts~logamazon_ts)
res_fit_ci<-fit_reg$residuals
CADFtest(res_fit_ci,type="drift", criterion = "BIC", max.lag.y=max.lag)

# The result of ADF test show that ADF(0)=-3.3749, and p-value=0.0123 which is 
# is smaller that the Engle-Granger ADF test statistics for
# one explanatory variable -3.14
# we reject H0 of no cointegration and conclude that are cointegrated.


##########################################################
#### Vector Error Correction  #####
logdata<-data.frame(log(nasdaq$IXIC.Close),log(amazon_data$AMZN.Close))
names(logdata)<-c("lognasdaq", "logamazon")
attach(logdata)

# Use Johansen's approach to test for cointegration
VARselect(logdata,lag.max =max.lag,type="const")
# According to the result of the Schwarz criterion selects
# the order 5 for the VAR on the time series in levels.
trace_test<-ca.jo(logdata,type="trace",K=5,ecdet="const",spec="transitory")
summary(trace_test)
maxeigen_test<-ca.jo(logdata,type="eigen",K=5,ecdet="const",spec="transitory")
summary(maxeigen_test)


# Estimate a Vector Error Correction Model
fit_vecm1<-cajorls(trace_test,r=1)
fit_vecm1

fit_vecm2<-cajorls(maxeigen_test,r=1)
fit_vecm2


# Use VECM to forecast the time series 
fit_var<-vec2var(trace_test,r=1)
myforecast<-predict(fit_var,n.ahead=365)
ts.plot(log(dfnasdaq))
lognasdaq_forecast<-ts(myforecast$fcst$lognasdaq[,1],frequency = 365,start=c(2018,12,19))
lognasdaq_lower<-ts(myforecast$fcst$lognasdaq[,2],frequency = 365,start=c(2018,12,19))
lognasdaq_upper<-ts(myforecast$fcst$lognasdaq[,3],frequency = 365, start=c(2018,12,19))
ts.plot(lognasdaq_forecast,lognasdaq_lower,lognasdaq_upper,col=c("black","red","red"))
title(main = "365-step-ahead forecast of log(nasdaq)")

plot.ts(logamazon)
logamazon_forecast<-ts(myforecast$fcst$logamazon[,1],frequency=365,start=c(2018,12))
logamazon_lower<-ts(myforecast$fcst$logamazon[,2],frequency=365,start=c(2018,12))
logamazon_upper<-ts(myforecast$fcst$logamazon[,3],frequency=365,start=c(2018,12))
ts.plot(logamazon_forecast,logamazon_lower,logamazon_upper,col=c("black","red","red"))
title(main = "365-step-ahead forecast of log(amazon)")


#########################################################################################









