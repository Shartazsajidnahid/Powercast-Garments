library(forecast)
setwd('/Users/nahid/Downloads/PowerGrid/code/PowerCast/code/tmp/')
X = read.csv('X.txt', header=F, stringsAsFactors=F)
X = as.numeric(X)
n_d_pred = read.csv('n_d_pred.txt', header=F, stringsAsFactors=F)
n_d_pred = as.numeric(n_d_pred)
u_freq = read.csv('u_freq.txt', header=F, stringsAsFactors=F)
u_freq = as.numeric(u_freq)
mod = auto.arima(ts(X, frequency=u_freq)) # for u vector
f = forecast(mod, h=n_d_pred) # for u vector - 2 days
write.table(f$mean, 'arima_Xp.txt', row.names=F, col.names=F)
