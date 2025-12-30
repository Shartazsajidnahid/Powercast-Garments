function [I_real_re, I_imag_re] = compforecast_fn(forecast_fn, I_real1, I_real2, I_imag1, I_imag2, ar_window, n_d_pred, u_freq)


if strcmp(forecast_fn, 'raw_AR')
    
    [I_real_re] = ar_fn2(ar_window,I_real1, I_real2);
    [I_imag_re] = ar_fn2(ar_window,I_imag1, I_imag2);
    
    
elseif strcmp(forecast_fn, 'raw_SAR')
    
    I_real_re = forecast_arima_seasonal_u2(I_real1, n_d_pred, u_freq);
    I_real_re = I_real_re';
    I_imag_re = forecast_arima_seasonal_u2(I_imag1, n_d_pred, u_freq);
    I_imag_re = I_imag_re';
    
    
end