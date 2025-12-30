function [M_new] = forecast_fn_tensor_SAR(M, n_d_pred)

r = size(M.U{1},2);
u_freq = 7;

for r_idx = 1: r
    pred_re_s = forecast_arima_seasonal_u2(M.U{1}(:,r_idx), n_d_pred, u_freq);
    M_U_pred(:,r_idx)  =  pred_re_s';
end


M_new.U{1} = M.U{1};
M_new.U{2} = M.U{2};
M_new.U{3} = M.U{3};
M_new.lambda = M.lambda;
M_new.U{1} = [M_new.U{1}; M_U_pred];

