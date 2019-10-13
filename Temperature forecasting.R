library(xts)
library(dygraphs)
library(zoo)
library(imputeTS)
library(forecast)
library(lubridate)
library(ggplot2)
library(scales)
library(prophet)
library(doParallel)
library(smooth)
library(aTSA)
library(MTS)
library(forecastHybrid)
library(stats)

cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

# Data preparation
getwd()
# read data and prepare for the analysis

getwd()
Analysis_starting_date = '2015-10-30 00:00:00'
k=7
Accuracy <- NULL

duom<-read.csv("../input/temperature.csv")
duom_new <- duom[,c(1,27)]
duom_new$Toronto <- duom_new$Toronto- 273.15
duom_new$datetime <- as.POSIXct(duom_new$datetime, '%Y-%m-%d %H:%M:%S', tz = "America/Toronto")
duom_new <- duom_new[duom_new$datetime >= Analysis_starting_date,]
duom_new <- duom_new[complete.cases(duom_new), ]

duom_wind <- read.csv("../input/humidity.csv")
duom_new_wind <- duom_wind[,c(1,27)]
duom_new_wind$datetime <- as.POSIXct(duom_new_wind$datetime, '%Y-%m-%d %H:%M:%S', tz = "America/Toronto")
duom_new_wind <- duom_new_wind[duom_new_wind$datetime >= Analysis_starting_date,]
duom_new_wind <- duom_new_wind[complete.cases(duom_new_wind), ]

ggplot(data = duom_new, aes(x = datetime, y = Toronto))+
  geom_line(color = "#00AFBB", size = 0.5) + ggtitle('Air temperature 2013-2017')+
  xlab('Date') + ylab('Temperature in Celsius')


ggplot(data = duom_new[duom_new$datetime >= '2016-11-30 12:00:00',], aes(x = datetime, y = Toronto))+
  geom_line(color = "#00AFBB", size = 0.5) + ggtitle('Last one year (2017) temperature graph')+
  xlab('Date') + ylab('Temperature in Celsius')


ggplot(data = duom_new[duom_new$datetime >= '2017-11-23 12:00:00',], aes(x = datetime, y = Toronto))+
  geom_line(color = "#00AFBB", size = 0.5) + ggtitle('Last one week temperature graph')+
  xlab('Date') + ylab('Temperature in Celsius')

plotNA.distribution(duom_new$Toronto)
statsNA(duom_new$Toronto)

outlier_values <- boxplot.stats(duom_new$Toronto)$out
outlier_values


proc_outliers <- (length(outlier_values)/length(duom_new$Toronto))
percent(proc_outliers)


D_1 <- duom_new$Toronto

summary(D_1)
summary<- summary(D_1)

IQR <- summary[5]- summary[2]

sprintf("Tarpkvantilinis plotis yra: %f",IQR)


Upper <- summary[5]+3*IQR
Lower <- summary[2]-3*IQR


D_1[D_1 > Upper] <- NA
D_1[D_1 < Lower] <- NA
plot(D_1, type="l")
plotNA.distribution(D_1)
statsNA(D_1)

monthplot(D_1, type="l")

if (length(D_1)<5000){
  shapiro.test(D_1)
}

h <- hist(D_1, breaks = 10, density = 10,
          col = "lightgray", xlab = "Temperature", main = "Histogram with Normal curve") 
xfit <- seq(min(D_1), max(D_1), length = 40) 
yfit <- dnorm(xfit, mean = mean(D_1), sd = sd(D_1)) 
yfit <- yfit * diff(h$mids[1:2]) * length(D_1) 

lines(xfit, yfit, col = "black", lwd = 2)

Acf(D_1)

Start = Sys.time()
es_accuracy <- NULL
date_dif <- (duom_new$datetime[nrow(duom_new)] - duom_new$datetime[nrow(duom_new)-24*14] )/k
train_start_date <- duom_new$datetime[nrow(duom_new)-24*14]
for (i in (1:k)){
  train <- duom_new[duom_new$datetime <= train_start_date,]
  test <- duom_new[duom_new$datetime >= train_start_date & duom_new$datetime <= (train_start_date+date_dif),]
  msts_power <- msts(train$Toronto, seasonal.periods = c(24,24*365.25),ts.frequency=48, start = decimal_date(as.POSIXct(Analysis_starting_date)))  
  
  holt <- forecast::forecast(HoltWinters(msts_power),h = nrow(test))
  accuracy(holt,test$Toronto)
  acc_mean <- accuracy(holt,test$Toronto)
  es_accuracy <- rbind(es_accuracy, HOLT_WINTERS=acc_mean[2,])
  train_start_date <- train_start_date + date_dif
  
  print(paste0("Finished fold:",i," in ",round(Sys.time()-Start,2)," seconds"))
}

es_avg_acc <- t(sapply(by(es_accuracy,rownames(es_accuracy),colMeans),identity))
Accuracy <- rbind(Accuracy,es_avg_acc)
print(paste0("Mean forecasting took:", round(Sys.time()-Start,2)," seconds"))

msts_power <- msts(train$Toronto, seasonal.periods = c(24,24*365.25),ts.frequency=88, start = decimal_date(as.POSIXct(Analysis_starting_date))) 

holt <- forecast::forecast(HoltWinters(msts_power),h = 5*nrow(test))

autoplot(holt)+geom_line(color = "#00AFBB", size = 0.5) + ggtitle('Temperature forecasting (exponential smoothing method)')+ xlab('Date') + ylab('Temperature in Celsius')


autoplot(holt)+geom_line(color = "#00AFBB", size = 0.5) + ggtitle('Temperature forecasting (mean method)')+ xlab('Date') + ylab('Temperature in Celsius') + xlim(end(holt$mean)-((end(holt$mean)-start(msts_power))/10),end(holt$mean))


acf(holt$residuals, na.action=na.exclude)


Start = Sys.time()
mean_accuracy <- NULL
date_dif <- (duom_new$datetime[nrow(duom_new)] - duom_new$datetime[nrow(duom_new)-24*14] )/k
train_start_date <- duom_new$datetime[nrow(duom_new)-24*14]
for (i in (1:k)){
  train <- duom_new[duom_new$datetime <= train_start_date,]
  test <- duom_new[duom_new$datetime >= train_start_date & duom_new$datetime <= (train_start_date+date_dif),]
  msts_power <- msts(train$Toronto, seasonal.periods = c(24,24*365.25), start = decimal_date(as.POSIXct(Analysis_starting_date)))
  mean_baseline <- meanf(msts_power,h = nrow(test))
  autoplot(mean_baseline)+geom_line(color = "#00AFBB", size = 0.5) + ggtitle('Temperature forecasting (mean method)')+ xlab('Date') + ylab('Temperature in Celsius')
  accuracy(mean_baseline,test$Toronto)
  acc_mean <- accuracy(mean_baseline,test$Toronto)
  mean_accuracy <- rbind(mean_accuracy, Mean_constant=acc_mean[2,])
  mean_accuracy
  train_start_date <- train_start_date + date_dif
  print(paste0("Finished fold:",i," in ",round(Sys.time()-Start,2)," seconds"))
}
mean_avg_acc <- t(sapply(by(mean_accuracy,rownames(mean_accuracy),colMeans),identity))
Accuracy <- rbind(Accuracy,mean_avg_acc)
print(paste0("Mean forecasting took:", round(Sys.time()-Start,2)," seconds"))

msts_power <- msts(train$Toronto, seasonal.periods = c(24, 24*365.25), start = decimal_date(as.POSIXct("2013-11-30 12:00:00")))

mean_baseline <- meanf(msts_power,h = nrow(test))

autoplot(mean_baseline)+geom_line(color = "#00AFBB", size = 0.5) + ggtitle('Temperature forecasting (mean method)')+ xlab('Date') + ylab('Temperature in Celsius')

autoplot(mean_baseline)+geom_line(color = "#00AFBB", size = 0.5) + ggtitle('Temperature forecasting (mean method)')+ xlab('Date') + ylab('Temperature in Celsius') + xlim(end(mean_baseline$mean)-((end(mean_baseline$mean)-start(msts_power))/10),end(mean_baseline$mean))

stationary.test(mean_baseline$residuals, method = "adf")

stationary.test(mean_baseline$residuals, method = "kpss")

acf(mean_baseline$residuals, na.action=na.exclude)

Box.test(mean_baseline$residuals,type="Ljung-Box")

Start = Sys.time()
naive_accuracy <- NULL
date_dif <- (duom_new$datetime[nrow(duom_new)] - duom_new$datetime[nrow(duom_new)-24*14] )/k
train_start_date <- duom_new$datetime[nrow(duom_new)-24*14]
for (i in (1:k)){
  train <- duom_new[duom_new$datetime <= train_start_date,]
  test <- duom_new[duom_new$datetime >= train_start_date & duom_new$datetime <= (train_start_date+date_dif),]
  msts_power <- msts(train$Toronto, seasonal.periods = c(24,24*365.25), start = decimal_date(as.POSIXct(Analysis_starting_date)))
  naive_baseline <- naive(msts_power,h = nrow(test))
  autoplot(naive_baseline)+geom_line(color = "#00AFBB", size = 0.5) + ggtitle('Temperature forecasting (naive method)')+ xlab('Date') + ylab('Temperature in Celsius')
  accuracy(naive_baseline,test$Toronto)
  acc_naive <- accuracy(naive_baseline,test$Toronto)
  naive_accuracy <- rbind(naive_accuracy, Naive=acc_naive[2,])
  naive_accuracy
  train_start_date <- train_start_date + date_dif
  print(paste0("Finished fold:",i," in ",round(Sys.time()-Start,2)," seconds"))
}
naive_avg_acc <- t(sapply(by(naive_accuracy,rownames(naive_accuracy),colMeans),identity))
Accuracy <- rbind(Accuracy,naive_avg_acc)
print(paste0("Mean forecasting took:", round(Sys.time()-Start,2)," seconds"))

msts_power <- msts(train$Toronto, seasonal.periods = c(24,169,24*365.25), start = decimal_date(as.POSIXct("2013-11-30 12:00:00")))
naive_baseline <- naive(msts_power,h = nrow(test))

autoplot(naive_baseline)+geom_line(color = "#00AFBB", size = 0.5) + ggtitle('Temperature forecasting (naive method)')+ xlab('Date') + ylab('Temperature in Celsius')

autoplot(naive_baseline)+geom_line(color = "#00AFBB", size = 0.5) + ggtitle('Temperature forecasting (mean method)')+ xlab('Date') + ylab('Temperature in Celsius') + xlim(end(naive_baseline$mean)-((end(naive_baseline$mean)-start(msts_power))/10),end(naive_baseline$mean))
stationary.test(naive_baseline$residuals, method = "adf")
stationary.test(naive_baseline$residuals, method = "kpss")
acf(naive_baseline$residuals, na.action=na.exclude)
Box.test(naive_baseline$residuals,type="Ljung-Box")



Start <- Sys.time()
Sys.time()
tbats_accuracy <- NULL
date_dif <- (duom_new$datetime[nrow(duom_new)] - duom_new$datetime[nrow(duom_new)-24*14] )/k
train_start_date <- duom_new$datetime[nrow(duom_new)-24*14]
for (i in (1:k)){
  train <- duom_new[duom_new$datetime <= train_start_date,]
  test <- duom_new[duom_new$datetime >= train_start_date & duom_new$datetime <= (train_start_date+date_dif),]
  msts_power <- msts(train$Toronto, seasonal.periods = c(24,24*365.25), start = decimal_date(as.POSIXct(Analysis_starting_date)))
  tbats_power <- forecast::tbats(msts_power)
  f_tbats <- forecast::forecast(tbats_power, h=nrow(test))
  autoplot(f_tbats) + ggtitle('Temperature forecasting (TBATS method)')+ xlab('Date') + ylab('Temperature in Celsius')
  acc_tbats <-  accuracy(f_tbats,test$Toronto)
  tbats_accuracy <- rbind(tbats_accuracy, TBats=acc_tbats[2,])
  tbats_accuracy
  train_start_date <- train_start_date + date_dif
  print(paste0("Finished fold:",i," in ",round(Sys.time()-Start,2)," seconds"))
}

tbats_avg_acc <- t(sapply(by(tbats_accuracy,rownames(tbats_accuracy),colMeans),identity))
Accuracy <- rbind(Accuracy, tbats_avg_acc)
print(paste0("TBATS forecasting took:", round(Sys.time()-Start,2)," seconds"))

msts_power <- msts(train$Toronto, seasonal.periods = c(24,24*365.25), start = decimal_date(as.POSIXct(Analysis_starting_date)))
tbats_power <- forecast::tbats(msts_power)
f_tbats <-  forecast::forecast(tbats_power, h = nrow(test))

autoplot(f_tbats) + ggtitle('Temperature forecasting (TBATS method)')+ xlab('Date') + ylab('Temperature in Celsius')
autoplot(f_tbats)+geom_line(color = "#00AFBB", size = 0.5) + ggtitle('Temperature forecasting (mean method)')+ xlab('Date') + ylab('Temperature in Celsius') + xlim(end(f_tbats$mean)-((end(f_tbats$mean)-start(msts_power))/10),end(f_tbats$mean))

stationary.test(as.numeric(tbats_power$errors), method = "adf")
stationary.test(as.numeric(tbats_power$errors), method = "kpss")
acf(tbats_power$errors, na.action=na.exclude)
Box.test(tbats_power$errors,type="Ljung-Box")




Start <- Sys.time()
Sys.time()
prophet_accuracy <- NULL
date_dif <- (duom_new$datetime[nrow(duom_new)] - duom_new$datetime[nrow(duom_new)-24*14] )/k
train_start_date <- duom_new$datetime[nrow(duom_new)-24*14]
for (i in (1:k)){
  train <- duom_new[duom_new$datetime <= train_start_date,]
  test <- duom_new[duom_new$datetime >= train_start_date & duom_new$datetime <= (train_start_date+date_dif),]
  colnames(train) <- c('ds','y')
  fit_prophet <- prophet(train, weekly.seasonality = FALSE, yearly.seasonality=TRUE)
  future_duq <- data.frame(test$datetime)
  colnames(future_duq) <- 'ds'
  f_prophet <- predict(fit_prophet,future_duq)
  plot(fit_prophet,f_prophet)
  prophet_plot_components(fit_prophet,f_prophet)
  
  acc_prophet <- accuracy(f_prophet$yhat, test$Toronto)
  
  prophet_accuracy <- rbind(prophet_accuracy, Prophet=acc_prophet[1,])
  
  train_start_date <- train_start_date + date_dif
  print(paste0("Finished fold:",i," in ",round(Sys.time()-Start,2)," seconds"))
}

prophet_avg_acc <- t(sapply(by(prophet_accuracy,rownames(prophet_accuracy),colMeans),identity))
prophet_avg_acc <- cbind(prophet_avg_acc, MASE=NA, ACF1=NA)
Accuracy <- rbind(Accuracy,prophet_avg_acc)
print(paste0("Mean forecasting took:", round(Sys.time()-Start,2)," seconds"))  


residuals <- f_prophet$yhat - test$Toronto
stationary.test(residuals, method = "adf")
stationary.test(residuals, method = "kpss")
acf(residuals, na.action=na.exclude)
Box.test(residuals,type="Ljung-Box")





Accuracy <- as.data.frame(Accuracy)
ggplot(Accuracy, aes(x=rownames(Accuracy), y=RMSE)) + xlab("Forecasting method")+
  geom_bar(stat="identity")+ggtitle("Average RMSE by methods with k fold cross validation")

ggplot(Accuracy , aes(x=rownames(Accuracy), y=MAE)) + xlab("Forecasting method")+
  geom_bar(stat="identity")+ggtitle("Average MAE by methods with k fold cross validation")






# FORECASTING ERROR BY HORIZON WIDTH

Start = Sys.time()
mean_accuracy_horizon <- NULL
mean_accur <- NULL

# Horizontu skaiicus
h = 14
# Stebejimu viename horizonte skaicius
n_in_h = 24

(duom_new$datetime[nrow(duom_new)] - duom_new$datetime[nrow(duom_new)-n_in_h*h])/h


date_dif <- (duom_new$datetime[nrow(duom_new)] - duom_new$datetime[nrow(duom_new)-n_in_h*h])/h

train <- duom_new[duom_new$datetime <= duom_new$datetime[nrow(duom_new)-n_in_h*h],]
train_prop <- train
train_wind <- duom_new_wind[duom_new_wind$datetime <= duom_new_wind$datetime[nrow(duom_new_wind)-n_in_h*h],]


for (i in (1:h) ){
  
  print("----------------------------------------")
  print(paste0("Starting horizon number:",i))
  print("----------------------------------------")
  test <- duom_new[duom_new$datetime >= duom_new$datetime[nrow(duom_new)-n_in_h*h] & duom_new$datetime <= (duom_new$datetime[nrow(duom_new)-n_in_h*h-1]+i*date_dif),]
  
  test_wind <- duom_new_wind[duom_new_wind$datetime >= duom_new_wind$datetime[nrow(duom_new_wind)-n_in_h*h] & duom_new_wind$datetime <= (duom_new_wind$datetime[nrow(duom_new_wind)-n_in_h*h-1]+i*date_dif),]
  # 
  # 
  # arimax_model <- auto.arima(train$Toronto, xreg=train_wind$Toronto)
  # summary(arimax_model)
  # f_fourier <-  forecast::forecast(arimax_model, xreg=test_wind$Toronto, h=nrow(test_wind))
  # acc_arimax <- accuracy(f_fourier,test$Toronto)
  # acc_arimax <- cbind(acc_arimax, i)
  # mean_accuracy_horizon <- rbind(mean_accuracy_horizon, ARIMAX=acc_arimax[2,])
  # 
  # 
  # 
  # duom_total_train <- merge(x =train, y = train_wind, by = "datetime", all = TRUE)
  # duom_total_test <- merge(x =test, y = test_wind, by = "datetime", all = TRUE)
  # duom_train <- duom_total_train[,2:3]
  # duom_test <- duom_total_test[,2:3]
  # duom_train <- as.matrix(sapply(duom_train, as.numeric)) 
  # varma_model <- MTS::VARMA(duom_train, p=2, q=2)
  # varma_pred  <- VARMApred(varma_model, h=nrow(test))
  # accuracy(as.ts(varma_pred$pred)[,1], duom_test$Toronto.x)
  # acc_varma <- accuracy(as.ts(varma_pred$pred)[,1], duom_test$Toronto.x)
  # acc_varma <- cbind(acc_varma, i)
  # mean_accuracy_horizon <- rbind(mean_accuracy_horizon, VARMA=acc_varma[1,])
  # 
  
  
  msts_power <- msts(train$Toronto, seasonal.periods = c(24,24*365.25), start = decimal_date(as.POSIXct("2013-11-30 12:00:00")))
  
  naive_baseline <- naive(msts_power,h = nrow(test))
  acc_naive <- accuracy(naive_baseline,test$Toronto)
  acc_naive <- cbind(acc_naive,i)
  mean_accuracy_horizon <- rbind(mean_accuracy_horizon, Naive=acc_naive[2,])
  print(paste0("Finished naive method. Time:",round(Sys.time()-Start,2))) 
  
  
  mean_baseline <- meanf(msts_power,h = nrow(test))
  acc_mean <- accuracy(mean_baseline,test$Toronto)
  acc_mean <- cbind(acc_mean,i)
  mean_accuracy_horizon <- rbind(mean_accuracy_horizon, Mean_constant=acc_mean[2,])
  print(paste0("Finished mean method. Time:",round(Sys.time()-Start,2))) 
  
  # 
  # snaive_baseline <- snaive(msts_power,h = nrow(test))
  # acc_snaive <- accuracy(snaive_baseline,test$Toronto)
  # acc_snaive <- cbind(acc_snaive,i)
  # mean_accuracy_horizon <- rbind(mean_accuracy_horizon, Snaive=acc_snaive[2,])
  # print(paste0("Finished snaive method. Time:",round(Sys.time()-Start,2))) 
  
  
  tbats_power <- forecast::tbats(msts_power)
  f_tbats <- forecast::forecast(tbats_power, h = nrow(test))
  acc_tbats <-  accuracy(f_tbats,test$Toronto)
  acc_tbats <- cbind(acc_tbats,i)
  mean_accuracy_horizon <- rbind(mean_accuracy_horizon, TBats=acc_tbats[2,])
  print(paste0("Finished TBATS method. Time:",round(Sys.time()-Start,2)))
  
  
  colnames(train_prop) <- c('ds','y')
  fit_prophet <- prophet(train_prop, weekly.seasonality = FALSE, daily.seasonality = FALSE)
  future_duq <- data.frame(test$datetime)
  colnames(future_duq) <- 'ds'
  f_prophet <- predict(fit_prophet,future_duq)
  acc_prophet <- accuracy(f_prophet$yhat, test$Toronto)
  acc_prophet <- cbind(acc_prophet, MASE=NA, ACF1=NA)
  acc_prophet <- cbind(acc_prophet,i)
  mean_accuracy_horizon <- rbind(mean_accuracy_horizon, Prophet_yearly=acc_prophet[1,])
  print(paste0("Finished Prophet (daily) method. Time:",round(Sys.time()-Start,2)))
  
  
  fit_prophet <- prophet(train_prop, weekly.seasonality = FALSE, yearly.seasonality = FALSE)
  future_duq <- data.frame(test$datetime)
  colnames(future_duq) <- 'ds'
  f_prophet <- predict(fit_prophet,future_duq)
  acc_prophet <- accuracy(f_prophet$yhat, test$Toronto)
  acc_prophet <- cbind(acc_prophet, MASE=NA, ACF1=NA)
  acc_prophet <- cbind(acc_prophet,i)
  mean_accuracy_horizon <- rbind(mean_accuracy_horizon, Prophet_yearly=acc_prophet[1,])
  print(paste0("Finished Prophet (daily) method. Time:",round(Sys.time()-Start,2)))
  
  
  fit_prophet <- prophet(train_prop, weekly.seasonality = FALSE, yearly.seasonality = TRUE)
  future_duq <- data.frame(test$datetime)
  colnames(future_duq) <- 'ds'
  f_prophet <- predict(fit_prophet,future_duq)
  acc_prophet <- accuracy(f_prophet$yhat, test$Toronto)
  acc_prophet <- cbind(acc_prophet, MASE=NA, ACF1=NA)
  acc_prophet <- cbind(acc_prophet,i)
  mean_accuracy_horizon <- rbind(mean_accuracy_horizon, Prophet_daily_yearly=acc_prophet[1,])
  print(paste0("Finished Prophet (daily) method. Time:",round(Sys.time()-Start,2)))
  
  print("----------------------------------------")
  print(paste0("Finished horizon:",i," in ",round(Sys.time()-Start,2)," seconds"))  
  print("----------------------------------------")
}

mean_accur<- as.data.frame(mean_accuracy_horizon)
mean_accur$Method <- rownames(mean_accur)


ggplot(mean_accur, aes(y=RMSE, x=i, group=Method, colour= Method))+geom_line(size = 0.5) + ggtitle('Forecasting error (RMSE) by forecasting horizon')+ xlab('Horizon') + ylab('RMSE')


ggplot(mean_accur, aes(y=MAE, x=i, group=Method, colour= Method))+geom_line(size = 0.5) + ggtitle('Forecasting error (RMSE) by forecasting horizon')+xlab('Horizon') + ylab('MAE')

print("----------------------------------------")
print(paste0("FINISHED ALL HORIZONS in", round(Sys.time()-Start,2)," seconds"))  
print("----------------------------------------")