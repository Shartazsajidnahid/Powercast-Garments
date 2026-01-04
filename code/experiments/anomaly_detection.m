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
dataname = 'Garments';

data_out_folder = '../../plots/';
plot_save_flag = 1;
fontsize_in  =15;
linewidth_in = 3;
%%
% % n_ds = 6;
n_ds = 5;
forecast_fns = {'tensor_AR'};

sigma = 0.5;
window_size = 5;
% % r = 2;
r = 1;
ar_window = 1;

for n_d_pred_idx = 1: length(n_ds)
    % n_d_pred = 2%1;
    n_d_pred = n_ds(n_d_pred_idx);
    
    
    
    if strcmp(dataname, 'CMU')
        n_d_removed = 4;
    elseif strcmp(dataname, 'Garments')
        n_d_removed = 3;
    elseif strcmp(dataname, 'LBNL')
        n_d_removed = 3;
    end
    
    % % % [V_real1, I_real1, V_imag1, I_imag1, myTensor, V_real2, I_real2, V_imag2, I_imag2, n_d,label_tmp_unit, label_tmp,label_tmp2, u_freq] = processData(...
    % % %     dataname, n_d_pred, n_d_removed, sigma, window_size, fontsize_in, plot_save_flag, linewidth_in,data_out_folder);
    
    
    [V_real1, I_real1, V_imag1, I_imag1, myTensor, V_real2, I_real2, V_imag2, I_imag2, n_d,label_tmp_unit, label_tmp,label_tmp2, u_freq] = processData(...
        dataname, n_d_pred, n_d_removed, sigma, window_size, fontsize_in, plot_save_flag, linewidth_in,data_out_folder);
    
    t_tot = length(V_real1);
    t_tot2 = length(V_real2);
    
    %% run parafac
    
    X = tensor(myTensor);
    M = parafac_als(X, r);
    % % % plot_parafac_res(M,n_d,label_tmp, fontsize_in, linewidth_in, plot_save_flag,data_out_folder,dataname);
    
    
    %% reconstruct tensor - check the fitting result
    re_myTensor = reconstruct_tensor(M, n_d, r);
    % plot_retensor_res(myTensor, re_myTensor);
    
    data_proc_method = 1;
    % let's restore the time sequences!
    % % [I_real_re_check, I_imag_re_check] = tensor2data(re_myTensor, V_real1, V_imag1, data_proc_method, dataname);
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
    % subplot(2,1,2);
    % plot(I_imag1,'k', 'linewidth',2); hold on;
    % plot(I_imag_re_check,':r', 'linewidth',2);
    % set(gca,'XTick',[1:24:t_tot]);
    % set(gca,'XTickLabel',label_tmp);
    % set(gca,'XTickLabelRotation',90);
    % for i = 1: ceil(t_tot/24)
    %     plot([(i-1)*24 (i-1)*24], [min(I_imag1) max(I_imag1)], ':k');
    % end
    
    err_check_re = computeErr(I_real1, I_real_re_check);
    % err_check_imag = computeErr(I_imag1, I_imag_re_check);
    
    
    
    
    
    %% forecast
    % elongate the M.U{1} vector and reconstruct the tensor
    
    
    
    for fns_idx = 1: length(forecast_fns)
        
        forecast_fn = forecast_fns{fns_idx};
        
        
        
        if strcmp(forecast_fn(1:3), 'ten')
            [I_real_re, I_imag_re, M_new] = myforecast_fn(forecast_fn, M, ar_window, n_d_pred, n_d, V_real2, V_imag2, dataname);
            % % plot_forecast_res(M_new,n_d,n_d_pred,I_real2,I_imag2,I_real_re,I_imag_re, label_tmp2, t_tot,forecast_fn, fontsize_in);
            plot_forecast_res(M_new,n_d,n_d_pred,I_real2,I_real2,I_real_re,I_real_re, label_tmp2, t_tot,forecast_fn, fontsize_in);
            
            err_check_re = computeErr(I_real2(t_tot+1:t_tot2), I_real_re(t_tot+1:t_tot2));
            % err_check_imag = computeErr(I_imag2(t_tot+1:t_tot2), I_imag_re(t_tot+1:t_tot2));
            
            
        elseif strcmp(forecast_fn(1:3), 'raw')
            [I_real_re, I_imag_re] = compforecast_fn(forecast_fn, I_real1, I_real2, I_imag1, I_imag2, ar_window, n_d_pred*u_freq, u_freq);
            nn = length(I_real1);
            plot_forecast_res2(nn, I_real_re, I_imag_re, I_real2, I_imag2,forecast_fn,label_tmp2, fontsize_in);
            
            err_check_re = computeErr(I_real2(t_tot+1:t_tot2), I_real_re);
            err_check_imag = computeErr(I_imag2(t_tot+1:t_tot2), I_imag_re);
            
        end
        if plot_save_flag
            saveas(gcf, [data_out_folder dataname '_' forecast_fn '_npred' num2str(n_d_pred) '_res'],'png');
        end
        
        
        
        %%% accumulate all the results over trials (fn, n_d_pred)
        err_all_re(fns_idx,n_d_pred_idx) = err_check_re;
        % err_all_imag(fns_idx,n_d_pred_idx) = err_check_imag;
        
        predictions_re{fns_idx,n_d_pred_idx} = I_real_re;
        % predictions_imag{fns_idx,n_d_pred_idx} = I_imag_re;
        
        
        
    end
    
    %     close all
end

err_re = abs(I_real2(t_tot+1:end) - I_real_re(t_tot+1:end)) ./ I_real2(t_tot+1:end);
% err_imag = abs(I_imag2(t_tot+1:end) - I_imag_re(t_tot+1:end)) ./ I_imag2(t_tot+1:end);

err_imag = 0; %% nahid


err_tot = err_re + err_imag;


cutoff_val = 0.2;

figure;
subplot(2,1,1); plot(err_re);
% % subplot(2,1,2); plot(err_imag);

figure;
plot(err_tot);hold on;
plot(1:length(err_tot), cutoff_val * ones(length(err_tot),1));

% err_timepts = find(err_tot >= cutoff_val);
err_timepts = find(err_tot <= cutoff_val);


figure;
subplot(2,1,1); hold on;
h1 = plot(I_real2(1:t_tot),'color',[0.5 0.5 0.5], 'linewidth',2); hold on;
h2 = plot((t_tot+1:t_tot2),I_real2(t_tot+1:end),'.k','markersize',5);
h3 = plot((t_tot+1:t_tot2),I_real_re(t_tot+1:end),'-r', 'linewidth',2);
set(gca,'XTick',[1:24:t_tot2]);
set(gca,'XTickLabel',label_tmp);
set(gca,'XTickLabelRotation',90);
for i = 1: ceil(t_tot2/24)
    plot([(i-1)*24 (i-1)*24], [min(I_real1) max(I_real1)], ':k');
end
f1 = [1 2 3 4];
v1 = [t_tot+1 min(I_real2); t_tot+1 max(I_real2) ;t_tot2 max(I_real2);t_tot2 min(I_real2)];
patch('Faces',f1, 'Vertices',v1, 'EdgeColor','none','FaceColor',[0 0 0.1],'FaceAlpha',.05);
for i = 1: length(err_timepts)
    f1 = [1 2 3 4];
    v1 = [err_timepts(i)+t_tot min(I_real1); err_timepts(i)+t_tot max(I_real1) ;err_timepts(i)+1+t_tot max(I_real1);err_timepts(i)+1+t_tot min(I_real1)];
    
    patch('Faces',f1, 'Vertices',v1, 'EdgeColor','none','FaceColor',[1 0 0],'FaceAlpha',.3);
    
end
ylabel('I_r [A]', 'fontsize', fontsize_in);
set(gca, 'XTickLabel', []);



% % % subplot(2,1,2);
% % % h1 = plot(I_imag2(1:t_tot),'color',[0.5 0.5 0.5],'linewidth',2); hold on;
% % % h2 = plot((t_tot+1:t_tot2),I_imag2(t_tot+1:end),'.k', 'markersize',5);
% % % h3 = plot((t_tot+1:t_tot2),I_imag_re(t_tot+1:end),'-r', 'linewidth',2);
% % % set(gca,'XTick',[1:24:t_tot2]);
% % % set(gca,'XTickLabel',label_tmp2);
% % % set(gca,'XTickLabelRotation',90);
% % % for i = 1: ceil(t_tot2/24)
% % %     plot([(i-1)*24 (i-1)*24], [min(I_imag1) max(I_imag1)], ':k');
% % % end
% % % f1 = [1 2 3 4];
% % % v1 = [t_tot+1 min(I_imag2); t_tot+1 max(I_imag2) ;t_tot2 max(I_imag2);t_tot2 min(I_imag2)];
% % % patch('Faces',f1, 'Vertices',v1, 'EdgeColor','none','FaceColor',[0 0 0.1],'FaceAlpha',.05);
% % % 
% % % for i = 1: length(err_timepts)
% % %     f1 = [1 2 3 4];
% % %     v1 = [err_timepts(i)+t_tot min(I_imag1); err_timepts(i)+t_tot max(I_imag1) ;err_timepts(i)+1+t_tot max(I_imag1);err_timepts(i)+t_tot+1 min(I_imag1)];
% % % 
% % %     patch('Faces',f1, 'Vertices',v1, 'EdgeColor','none','FaceColor',[1 0 0],'FaceAlpha',.3);
% % % 
% % % end

legend([h2 h3],{'Truth','PowerCast'},'orientation','horizontal','location','none', 'position', [0.10826 0.46839 0.81217 0.06990], 'fontsize',fontsize_in);
ylabel('I_i [A]', 'fontsize', fontsize_in);

if 1
    saveas(gcf, [data_out_folder dataname '_anomaly2'],'png');
end