function myTensor = construct_tensor(V_real,I_real,V_imag,I_imag,n_d,t_d,data_proc_method,window_size,sigma, dataname)


% myTensor = zeros(n_d,t_d,4);
myTensor = zeros(n_d,t_d,2);

% data_proc_method = 1;
% window_size = 5; % should be odd number
% sigma = 0.5;


for day_idx = 1: n_d
    V_real_in = V_real((day_idx-1)*t_d+1:(day_idx-1)*t_d+t_d);
    I_real_in = I_real((day_idx-1)*t_d+1:(day_idx-1)*t_d+t_d);
    % V_imag_in = V_imag((day_idx-1)*t_d+1:(day_idx-1)*t_d+t_d);
    I_imag_in = I_imag((day_idx-1)*t_d+1:(day_idx-1)*t_d+t_d);
    
    if data_proc_method == 1 % CMU data
        for tt = 1: t_d
            
            V_real_in_tmp = windowing_fn(V_real_in, tt, t_d, window_size, sigma);
            % V_imag_in_tmp = windowing_fn(V_imag_in, tt, t_d, window_size, sigma);
            I_real_in_tmp = windowing_fn(I_real_in, tt, t_d, window_size, sigma);
            I_imag_in_tmp = windowing_fn(I_imag_in, tt, t_d, window_size, sigma);
            
            
            if strcmp(dataname, 'CMU')
                [b, ~, ~] = LinFit(V_real_in_tmp, I_real_in_tmp);
                fit_alpha_r = b(1);
                fit_G = b(2);
                [b, ~, ~] = LinFit(V_real_in_tmp, I_imag_in_tmp);
                fit_alpha_i = b(1);
                fit_B = b(2);
            elseif strcmp(dataname, 'Garments')
                [b, ~, ~] = LinFit(V_real_in_tmp, I_real_in_tmp);
                if isnan(b(1))
                    fit_alpha_r = 0;
                else
                    fit_alpha_r = b(1);
                end

                if isnan(b(2))
                    fit_G = 0;
                else
                    fit_G = b(2);
                end

                [b, ~, ~] = LinFit(V_real_in_tmp, I_imag_in_tmp);
                if isnan(b(1))
                    fit_alpha_i = 0;
                else
                    fit_alpha_i = b(1);
                end

                if isnan(b(2))
                    fit_B = 0;
                else
                    fit_B = b(2);
                end
                
                % [b, ~, ~] = LinFit(V_real_in_tmp, I_imag_in_tmp);
                % fit_alpha_i = b(1);
                % fit_B = b(2);
            elseif strcmp(dataname, 'LBNL')
                bic_param = 2^3;
                [b, ~, ~, ~] = BIGFitSegment(V_real_in_tmp, V_imag_in_tmp, I_real_in_tmp, I_imag_in_tmp, 0, bic_param);
                fit_G = b(1);
                fit_B = b(2);
                fit_alpha_r = b(3);
                fit_alpha_i = b(4);
            end
            
            myTensor(day_idx, tt,1) = fit_G;
            myTensor(day_idx, tt,2) = fit_alpha_r;
            myTensor(day_idx, tt,3) = fit_B;
            myTensor(day_idx, tt,4) = fit_alpha_i;
        end
        
    elseif data_proc_method == 2
        seg_method = 'greedy';
        [fit_G, fit_alpha_r, fit_B, fit_alpha_i, yhat_cut_auto_f, I_real_hat, I_imag_hat] = run_alg(seg_method, V_real_in, V_imag_in, I_real_in, I_imag_in);
        
        for seg_idx = 1:size(yhat_cut_auto_f,2)
            myTensor(day_idx, yhat_cut_auto_f{seg_idx},1) = fit_G(seg_idx);
            myTensor(day_idx, yhat_cut_auto_f{seg_idx},2) = fit_alpha_r(seg_idx);
            myTensor(day_idx, yhat_cut_auto_f{seg_idx},3) = fit_B(seg_idx);
            myTensor(day_idx, yhat_cut_auto_f{seg_idx},4) = fit_alpha_i(seg_idx);
        end
        
        
    end
    
    
    close all
end

