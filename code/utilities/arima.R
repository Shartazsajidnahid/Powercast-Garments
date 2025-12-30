library(forecast)
setwd('~/Desktop/bryan-papers/heart/code')
X = read.csv('../temp/X.txt', header=F, stringsAsFactors=F)$V1
X = as.numeric(X)
mod = auto.arima(ts(X))
f = forecast(mod, h=1000)
write.table(f$mean, '../temp/arima_Xp.txt', row.names=F, col.names=F)
