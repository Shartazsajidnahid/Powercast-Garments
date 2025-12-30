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

data_out_folder = '../../plots/';
plot_save_flag = 1;
fontsize_in =15;
linewidth_in = 3;


%%
% number of days to forecast
n_ds = 1:2;

% methods
forecast_fns = {'tensor_AR', 'tensor_SAR', 'raw_AR', 'raw_SAR'};

% Gaussian parameters
sigma = 0.5;
window_size = 5;
% rank
r = 2;
% AR(p)
ar_window = 1;



for n_d_pred_idx = 1: length(n_ds)
    n_d_pred = n_ds(n_d_pred_idx);
    
    
    if strcmp(dataname, 'CMU')
        n_d_removed = 4;
    elseif strcmp(dataname, 'LBNL')
        n_d_removed = 3;
    end
    
    [V_real1, I_real1, V_imag1, I_imag1, myTensor, V_real2, I_real2, V_imag2, I_imag2, n_d,label_tmp_unit, label_tmp,label_tmp2, u_freq] = processData(...
        dataname, n_d_pred, n_d_removed, sigma, window_size, fontsize_in, plot_save_flag, linewidth_in,data_out_folder);
    
    t_tot = length(V_real1);
    t_tot2 = length(V_real2);
    
    %% run parafac
    X = tensor(myTensor);
    M = parafac_als(X, r);
    plot_parafac_res(M,n_d,label_tmp, fontsize_in, linewidth_in, plot_save_flag,data_out_folder,dataname);
    
    
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
        
        if strcmp(forecast_fn(1:3), 'ten')
            [I_real_re, I_imag_re, M_new] = myforecast_fn(forecast_fn, M, ar_window, n_d_pred, n_d, V_real2, V_imag2, dataname);
            plot_forecast_res(M_new,n_d,n_d_pred,I_real2,I_imag2,I_real_re,I_imag_re, label_tmp2, t_tot,forecast_fn, fontsize_in);
            
            err_check_re = computeErr(I_real2(t_tot+1:t_tot2), I_real_re(t_tot+1:t_tot2));
            err_check_imag = computeErr(I_imag2(t_tot+1:t_tot2), I_imag_re(t_tot+1:t_tot2));
            
            
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
        err_all_imag(fns_idx,n_d_pred_idx) = err_check_imag;
        
        predictions_re{fns_idx,n_d_pred_idx} = I_real_re;
        predictions_imag{fns_idx,n_d_pred_idx} = I_imag_re;
        
        
        
    end
    
    close all
end




for i = 1: n_d_pred
    x_label1{i} =[ num2str(i) ' day(s)'];
end
x_label = categorical(x_label1);

error_all_sum = err_all_re + err_all_imag;
figure;
bar(x_label, error_all_sum');
% legend({'Tensor-AR','Tensor-SAR','Raw-AR','Raw-SAR'}, 'location', 'eastoutside', 'fontsize',fontsize_in);
legend({'PowerCast','PowerCast-S','AR','AR-S'}, 'location', 'eastoutside', 'fontsize',fontsize_in);
xlabel('number of days for forecast', 'fontsize',fontsize_in);
ylabel('error sum of I_r and I_i', 'fontsize',fontsize_in);

if plot_save_flag
    saveas(gcf, [data_out_folder dataname '_err_per_days_res'],'png');
end


% individual for each day
hFig = figure;
p1X = [1]; p1Y = error_all_sum(1,1);
p2X = [2]; p2Y = error_all_sum(2,1);
p3X = [3]; p3Y = error_all_sum(3,1);
p4X = [4]; p4Y = error_all_sum(4,1);
p1 = bar(p1X,p1Y);
hold on;
p2 = bar(p2X,p2Y);
hold on;
p3 = bar(p3X,p3Y);
hold on;
p4 = bar(p4X,p4Y);
set(p1,'FaceColor','r');
set(p2,'FaceColor','m');
set(p3,'FaceColor',[0 0 0.5]);
set(p4,'FaceColor',[0 0.5 0]);
set(gca,'xticklabel',{})
set(gca, 'XTick', [1:4]);
xin{1} = 'PowerCast';
xin{2} = 'PowerCast-S';
xin{3} = 'AR';
xin{4} = 'SAR';
set(gca, 'XTickLabel', xin);
set(gca,'XTickLabelRotation',60);
set(gca, 'fontsize', fontsize_in+20);

ylabel('error sum of I_r and I_i', 'fontsize',fontsize_in+20);
set(hFig, 'Position', [100 100 500 800])
if plot_save_flag
    saveas(gcf, [data_out_folder dataname '_err_per_days_res_1dayonly'],'png');
end

improvement1 = (min(error_all_sum(4,1),error_all_sum(3,1)) - min(error_all_sum(2,1), error_all_sum(1,1))) / min(error_all_sum(4,1),error_all_sum(3,1));



% individual for each day
hFig = figure;
p1X = [1]; p1Y = error_all_sum(1,2);
p2X = [2]; p2Y = error_all_sum(2,2);
p3X = [3]; p3Y = error_all_sum(3,2);
p4X = [4]; p4Y = error_all_sum(4,2);
p1 = bar(p1X,p1Y);
hold on;
p2 = bar(p2X,p2Y);
hold on;
p3 = bar(p3X,p3Y);
hold on;
p4 = bar(p4X,p4Y);
set(p1,'FaceColor','red');
set(p2,'FaceColor','m');
set(p3,'FaceColor',[0 0 0.5]);
set(p4,'FaceColor',[0 0.5 0]);
set(gca,'xticklabel',{})
set(gca, 'XTick', [1:4]);
xin{1} = 'PowerCast';
xin{2} = 'PowerCast-S';
xin{3} = 'AR';
xin{4} = 'SAR';
set(gca, 'XTickLabel', xin);
set(gca,'XTickLabelRotation',60);
set(gca, 'fontsize', fontsize_in+20);

ylabel('error sum of I_r and I_i', 'fontsize',fontsize_in+20);
set(hFig, 'Position', [100 100 500 800])
if plot_save_flag
    saveas(gcf, [data_out_folder dataname '_err_per_days_res_2dayonly'],'png');
end

improvement2 = (min([error_all_sum(4,2),error_all_sum(3,2)]) - min([error_all_sum(2,2), error_all_sum(1,2)])) / min([error_all_sum(4,2),error_all_sum(3,2)]);



plot_forecast_res_all(n_ds, predictions_re, predictions_imag, I_real2, I_imag2, label_tmp2, forecast_fns, dataname, data_out_folder,plot_save_flag, fontsize_in);


