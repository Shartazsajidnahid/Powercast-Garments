function myTensor = construct_tensor_segmentation(V_real,I_real,V_imag,I_imag,n_d,t_d,method_funcs,method_idx, BIC_params, trial_idx)


% myTensor = zeros(n_d,t_d,4);



% data_proc_method = 1;
% window_size = 5; % should be odd number
% sigma = 0.5;


for day_idx = 1: n_d
    V_real_in(day_idx,:) = V_real((day_idx-1)*t_d+1:(day_idx-1)*t_d+t_d);
    I_real_in(day_idx,:) = I_real((day_idx-1)*t_d+1:(day_idx-1)*t_d+t_d);
    V_imag_in(day_idx,:) = V_imag((day_idx-1)*t_d+1:(day_idx-1)*t_d+t_d);
    I_imag_in(day_idx,:) = I_imag((day_idx-1)*t_d+1:(day_idx-1)*t_d+t_d);
end


[b_all, Irhat, Iihat, segs] = tensor_LinFitGreedy(V_real_in, V_imag_in, I_real_in, I_imag_in, method_funcs{method_idx}, BIC_params(trial_idx));

nseg = size(segs,2);

myTensor = zeros(n_d,nseg,4);
for j = 1: n_d
    for i = 1: nseg
        myTensor(j,i,1)  = b_all{j,i}(1); % G
        myTensor(j,i,3)  = b_all{j,i}(2); % B
        myTensor(j,i,2)  = b_all{j,i}(3); % ar
        myTensor(j,i,4)  = b_all{j,i}(4); % ai
    end
end

% % %     if data_proc_method == 1 % CMU data
% % %         for tt = 1: t_d
% % %
% % %             V_real_in_tmp = windowing_fn(V_real_in, tt, t_d, window_size, sigma);
% % %             V_imag_in_tmp = windowing_fn(V_imag_in, tt, t_d, window_size, sigma);
% % %             I_real_in_tmp = windowing_fn(I_real_in, tt, t_d, window_size, sigma);
% % %             I_imag_in_tmp = windowing_fn(I_imag_in, tt, t_d, window_size, sigma);
% % %
% % %
% % %             if strcmp(dataname, 'CMU')
% % %                 [b, ~, ~] = LinFit(V_real_in_tmp, I_real_in_tmp);
% % %                 fit_alpha_r = b(1);
% % %                 fit_G = b(2);
% % %                 [b, ~, ~] = LinFit(V_real_in_tmp, I_imag_in_tmp);
% % %                 fit_alpha_i = b(1);
% % %                 fit_B = b(2);
% % %             elseif strcmp(dataname, 'LBNL')
% % %                 bic_param = 2^3;
% % %                 [b, ~, ~, ~] = BIGFitSegment(V_real_in_tmp, V_imag_in_tmp, I_real_in_tmp, I_imag_in_tmp, 0, bic_param);
% % %                 fit_G = b(1);
% % %                 fit_B = b(2);
% % %                 fit_alpha_r = b(3);
% % %                 fit_alpha_i = b(4);
% % %             end
% % %
% % %             myTensor(day_idx, tt,1) = fit_G;
% % %             myTensor(day_idx, tt,2) = fit_alpha_r;
% % %             myTensor(day_idx, tt,3) = fit_B;
% % %             myTensor(day_idx, tt,4) = fit_alpha_i;
% % %         end
% % %
% % %     elseif data_proc_method == 2
% % %         seg_method = 'greedy';
% % %         [fit_G, fit_alpha_r, fit_B, fit_alpha_i, yhat_cut_auto_f, I_real_hat, I_imag_hat] = run_alg(seg_method, V_real_in, V_imag_in, I_real_in, I_imag_in);
% % %
% % %         for seg_idx = 1:size(yhat_cut_auto_f,2)
% % %             myTensor(day_idx, yhat_cut_auto_f{seg_idx},1) = fit_G(seg_idx);
% % %             myTensor(day_idx, yhat_cut_auto_f{seg_idx},2) = fit_alpha_r(seg_idx);
% % %             myTensor(day_idx, yhat_cut_auto_f{seg_idx},3) = fit_B(seg_idx);
% % %             myTensor(day_idx, yhat_cut_auto_f{seg_idx},4) = fit_alpha_i(seg_idx);
% % %         end
% % %
% % %
% % %     end
% % %
% % %
% % %     close all
end





