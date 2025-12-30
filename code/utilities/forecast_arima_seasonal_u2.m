function Xp = forecast_arima_seasonal_u2(X, n_d_pred,u_freq)

dlmwrite('../tmp/X.txt',X');
dlmwrite('../tmp/n_d_pred.txt',n_d_pred);
dlmwrite('../tmp/u_freq.txt',u_freq);

system('/usr/local/bin/R CMD BATCH arima_seasonal_u2.R');
Xp = dlmread('../tmp/arima_Xp.txt')';
% assert(length(Xp) == pred_len);
