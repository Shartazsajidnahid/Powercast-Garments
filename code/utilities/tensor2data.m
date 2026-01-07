function [I_real_re, I_imag_re] = tensor2data(myTensor, V_real, V_imag, data_proc_method, dataname)

% % % I_imag_re = 1; %% nahid

n_d = size(myTensor,1);
t_d = size(myTensor,2);

k = 1;
for day_idx = 1: n_d
    
    if data_proc_method == 1
        for tt = 1: t_d
            % retrieve parameters
            G_now = myTensor(day_idx, tt,1);
            alpha_r_now = myTensor(day_idx, tt,2);
            B_now = myTensor(day_idx, tt,3);
            alpha_i_now = myTensor(day_idx, tt,4);
            
            if strcmp(dataname, 'CMU')
                % recover I
                I_real_re(k,1) = alpha_r_now + G_now * V_real(k,1);
                I_imag_re(k,1) = alpha_i_now + B_now * V_real(k,1);
                k = k + 1;
            elseif strcmp(dataname, 'Garments')
                % recover I
                I_real_re(k,1) = alpha_r_now + G_now * V_real(k,1);
                I_imag_re(k,1) = alpha_i_now + B_now * V_real(k,1);
                k = k + 1;
            elseif strcmp(dataname, 'LBNL')
                I_real_re(k,1) = alpha_r_now + G_now * V_real(k,1) - B_now * V_imag(k,1);
                I_imag_re(k,1) = alpha_i_now + B_now * V_real(k,1) + G_now * V_imag(k,1);
                k = k + 1;
            end
        end
        
        
    end
end

