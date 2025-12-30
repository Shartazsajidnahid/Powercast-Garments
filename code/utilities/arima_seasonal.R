library(forecast)
setwd('/Users/hyunah/SVN_repo/hyunahs-code/PowerGrid/code/')
#X = read.csv('X.txt', header=F, stringsAsFactors=F)$V1
X = read.csv('X.txt', header=F, stringsAsFactors=F)
X = as.numeric(X)
#mod = auto.arima(ts(X, frequency=7)) # for u vector
mod = auto.arima(ts(X, frequency=24)) # for raw data
#mod = auto.arima(ts(X),seasonal=TRUE)
#f = forecast(mod, h=1000)
f = forecast(mod, h=48) # for raw data - 2 days * 24 hrs
#f = forecast(mod, h=2) # for u vector - 2 days
#write.table(f$mean, 'arima_Xp.txt', row.names=F, col.names=F)
write.table(f$mean, 'arima_Xp.txt', row.names=F, col.names=F)
