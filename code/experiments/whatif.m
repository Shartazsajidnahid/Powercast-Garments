clear all
close all
clc

rngseed = 1;
rng(rngseed);

includedir = '../../tensor_toolbox/';
addpath(genpath(includedir));
addpath(genpath('../utilities/'));
% choose data
dataname = 'CMU';
% dataname = 'LBNL';

data_out_folder = '../../plots/';
plot_save_flag = 1;
fontsize_in = 15;
linewidth_in = 3;
%%
n_ds = 2;
forecast_fns = {'tensor_AR'};

sigma = 0.5;
window_size = 5;
r = 2;
 ar_window = 1;


for n_d_pred_idx = 1: length(n_ds)
    
    n_d_pred = n_ds(n_d_pred_idx);
    
    
    
    if strcmp(dataname, 'CMU')
        n_d_removed = 4;
    elseif strcmp(dataname, 'LBNL')
        n_d_removed = 3;
    end
    
    [V_real1, I_real1, V_imag1, I_imag1, myTensor, V_real2, I_real2, V_imag2, I_imag2, n_d,label_tmp_unit, label_tmp,label_tmp2, u_freq] = processData(...
        dataname, n_d_pred, n_d_removed, sigma, window_size, fontsize_in, 0, linewidth_in,data_out_folder);
    
    t_tot = length(V_real1);
    t_tot2 = length(V_real2);
    
    %% run parafac
    
    X = tensor(myTensor);
    M = parafac_als(X, r);
    plot_parafac_res(M,n_d,label_tmp,fontsize_in,linewidth_in, 0,data_out_folder,dataname)
    
    
    
    %% reconstruct tensor - check the fitting result
    re_myTensor = reconstruct_tensor(M, n_d, r);
    plot_retensor_res(myTensor, re_myTensor);
    
    data_proc_method = 1;
    % let's restore the time sequences!
    [I_real_re_check, I_imag_re_check] = tensor2data(re_myTensor, V_real1, V_imag1, data_proc_method, dataname);
    
    % plot
    figure;
    subplot(2,1,1);
    plot(I_real1,'k', 'linewidth',2); hold on;
    plot(I_real_re_check,':r', 'linewidth',2);
    set(gca,'XTick',[1:24:t_tot]);
    set(gca,'XTickLabel',label_tmp);
    set(gca,'XTickLabelRotation',90);
    for i = 1: ceil(t_tot/24)
        plot([(i-1)*24 (i-1)*24], [min(I_real1) max(I_real1)], ':k');
    end
    subplot(2,1,2);
    plot(I_imag1,'k', 'linewidth',2); hold on;
    plot(I_imag_re_check,':r', 'linewidth',2);
    set(gca,'XTick',[1:24:t_tot]);
    set(gca,'XTickLabel',label_tmp);
    set(gca,'XTickLabelRotation',90);
    for i = 1: ceil(t_tot/24)
        plot([(i-1)*24 (i-1)*24], [min(I_imag1) max(I_imag1)], ':k');
    end
    
    err_check_re = computeErr(I_real1, I_real_re_check);
    err_check_imag = computeErr(I_imag1, I_imag_re_check);
    
    
    
    
    
    %% forecast
    % elongate the M.U{1} vector and reconstruct the tensor
   
    
    
    for fns_idx = 1: length(forecast_fns)
        
        forecast_fn = forecast_fns{fns_idx};
        
        
        
        % % %         if strcmp(forecast_fn(1:3), 'ten')
        [I_real_re, I_imag_re, M_new, re_myTensor_prd] = myforecast_fn(forecast_fn, M, ar_window, n_d_pred, n_d, V_real2, V_imag2, dataname);
        plot_forecast_res(M_new,n_d,n_d_pred,I_real2,I_imag2,I_real_re,I_imag_re, label_tmp2, t_tot,forecast_fn, fontsize_in);
        
        % what-if scenario simulation
        
        if strcmp(dataname, 'CMU')
            
            try_sets{1,1} = [1.1 1.2 1.3];
            try_sets{1,2} = 1;
            
            try_sets{2,1} = 1;
            try_sets{2,2} = [1.1 1.2 1.3];
            
        elseif strcmp(dataname, 'LBNL')
           
            
            try_sets{1,1} = 1;
            try_sets{1,2} = [3 5 10];
        end
        
        for tt = 1: size(try_sets,1)
            
            Gfactors = try_sets{tt,1};
            Bfactors = try_sets{tt,2};
            
            for g_idx = 1: length(Gfactors)
                for b_idx = 1: length(Bfactors)
                    Gfactor = Gfactors(g_idx);
                    Bfactor = Bfactors(b_idx);
                    [I_real_re_whatif{g_idx, b_idx}, I_imag_re_whatif{g_idx, b_idx}] = tensor2data_whatif(re_myTensor_prd, V_real2, V_imag2, data_proc_method, dataname, Gfactor, Bfactor);
                    
                    [vvv, iii] = max(abs(I_real2(t_tot:end) - I_real_re_whatif{g_idx, b_idx}(t_tot:end)));
                    times_more_real{tt}(g_idx, b_idx) =  vvv / I_real2(t_tot+iii-1);
                    
                    [vvv, iii] = max(abs(I_imag2(t_tot:end) - I_imag_re_whatif{g_idx, b_idx}(t_tot:end)));
                    times_more_imag{tt}(g_idx, b_idx) =  vvv / I_imag2(t_tot+iii-1);
                end
            end
            plot_forecast_whatif_res(I_real2,I_imag2,I_real_re,I_imag_re, I_real_re_whatif,I_imag_re_whatif, label_tmp2, t_tot,forecast_fn,Gfactors,Bfactors, fontsize_in);
            
            if plot_save_flag
                saveas(gcf, [data_out_folder dataname '_whatif_d' num2str(n_d_pred) '_' forecast_fn '_t' num2str(tt)],'png');
            end
            
            plot_forecast_whatif_mag_res(I_real2,I_imag2,I_real_re,I_imag_re, I_real_re_whatif,I_imag_re_whatif, label_tmp2, t_tot,forecast_fn,Gfactors,Bfactors, fontsize_in, plot_save_flag,data_out_folder, dataname,n_d_pred,tt);
            
        end
        
        
    end
    
    %     close all
end


