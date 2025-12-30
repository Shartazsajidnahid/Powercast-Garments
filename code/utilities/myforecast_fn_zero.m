function [I_real_re, I_imag_re, M_new] = myforecast_fn_zero(forecast_fn, M, ar_window, n_d_pred, n_d, V_real2, V_imag2, dataname)


r = size(M.U{1},2);

if strcmp(forecast_fn, 'tensor_AR')
    
    [M_new] = forecast_fn_tensor_AR(M, ar_window, n_d_pred, n_d);
    
    
elseif strcmp(forecast_fn, 'tensor_SAR')
    
    [M_new] = forecast_fn_tensor_SAR(M, n_d_pred);
    
end

% reconstruct tensor
re_myTensor_prd_tmp = reconstruct_tensor(M_new, n_d+n_d_pred, r);


re_myTensor_prd(:,:,1) = re_myTensor_prd_tmp(:,:,1);
re_myTensor_prd(:,:,2) = zeros(size(re_myTensor_prd_tmp(:,:,1)));
re_myTensor_prd(:,:,3) = re_myTensor_prd_tmp(:,:,2);
re_myTensor_prd(:,:,4) = zeros(size(re_myTensor_prd_tmp(:,:,1)));



% recover the I_real and I_imag from V_real and re_myTensor_prd
data_proc_method = 1;
[I_real_re, I_imag_re] = tensor2data(re_myTensor_prd, V_real2, V_imag2, data_proc_method, dataname);










