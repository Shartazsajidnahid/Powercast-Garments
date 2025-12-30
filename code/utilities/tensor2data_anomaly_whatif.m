function [I_real_re, I_imag_re] = tensor2data_anomaly_whatif(myTensor, V_real, V_imag, dataname,err_timepts,times_factor, t_tot)


n_d = size(myTensor,1);
t_d = size(myTensor,2);

G_now = reshape(transpose(myTensor(:,:,1)), n_d*t_d,1);
alpha_r_now = reshape(transpose(myTensor(:,:,2)), n_d*t_d,1);
B_now = reshape(transpose(myTensor(:,:,3)), n_d*t_d,1);
alpha_i_now = reshape(transpose(myTensor(:,:,4)), n_d*t_d,1);

if strcmp(dataname, 'CMU')
    % recover I
    I_real_re = alpha_r_now + G_now .* V_real;
    I_imag_re = alpha_i_now + B_now .* V_real;
elseif strcmp(dataname, 'LBNL')
    I_real_re = alpha_r_now + G_now .* V_real - B_now * V_imag;
    I_imag_re = alpha_i_now + B_now .* V_real + G_now * V_imag;
end

if strcmp(dataname, 'CMU')
    for i = 1: length(err_timepts)
         timenow = t_tot + err_timepts(i);
        I_real_re(timenow) = alpha_r_now(timenow) + times_factor * G_now(timenow) .* V_real(timenow);
        I_imag_re(timenow) = alpha_i_now(timenow) + times_factor * B_now(timenow).* V_real(timenow);
    end
elseif strcmp(dataname, 'LBNL')
    for i = 1: length(err_timepts)
        timenow = t_tot + err_timepts(i);
        I_real_re(timenow)  = alpha_r_now(timenow) + times_factor * G_now(timenow)  .* V_real(timenow)  - times_factor * B_now(timenow)  * V_imag(timenow) ;
        I_imag_re(timenow)  = alpha_i_now(timenow) + times_factor * B_now(timenow)  .* V_real(timenow)  + times_factor * G_now(timenow) * V_imag(timenow) ;
    end
end
