clear all
close all
clc

rngseed = 1;
rng(rngseed);

includedir = '../../tensor_toolbox/';

addpath(genpath(includedir));
addpath(genpath('../utilities/'));

dataname = 'CMU';

data_out_folder = '../../plots/';
plot_save_flag = 0;
fontsize_in  =15;
linewidth_in = 3;
%%
n_ds = 1;

forecast_fn = {'tensor_AR'};

sigma = 0.5;
window_size = 5;

n_d_pred = 1;

n_d_removeds = 18:-2:1;

for n_d_removed_idx = 1:length(n_d_removeds)
    
    n_d_removed = n_d_removeds(n_d_removed_idx);
    
    
    [V_real1, I_real1, V_imag1, I_imag1, myTensor, V_real2, I_real2, V_imag2, I_imag2, n_d,label_tmp_unit, label_tmp,label_tmp2, u_freq] = processData(...
        dataname, n_d_pred, n_d_removed, sigma, window_size, fontsize_in, plot_save_flag, linewidth_in,data_out_folder);
    close all
    
    t_tot = length(V_real1);
    t_tot2 = length(V_real2);

    tic;
    %% run parafac
    r = 2;
    X = tensor(myTensor);
    M = parafac_als(X, r);
    
    %% reconstruct tensor - check the fitting result
    re_myTensor = reconstruct_tensor(M, n_d, r);
    
    data_proc_method = 1;
    % let's restore the time sequences!
    [I_real_re_check, I_imag_re_check] = tensor2data(re_myTensor, V_real1, V_imag1, data_proc_method, dataname);
    
    
    %% forecast
    % elongate the M.U{1} vector and reconstruct the tensor
    ar_window = 1;
    
    [I_real_re, I_imag_re, M_new] = myforecast_fn(forecast_fn, M, ar_window, n_d_pred, n_d, V_real2, V_imag2, dataname);
    
    err_check_re = computeErr(I_real2(t_tot+1:t_tot2), I_real_re(t_tot+1:t_tot2));
    err_check_imag = computeErr(I_imag2(t_tot+1:t_tot2), I_imag_re(t_tot+1:t_tot2));
    
    times(n_d_removed_idx) = toc;
    
    n_pts(n_d_removed_idx) = t_tot;
    
   
end

save scalability_res n_pts times

load scalability_res n_pts times

figure; hold on;
plot(n_pts(2:end-1), times(2:end-1), '.k','markersize', 20, 'linewidth',2);
plot([n_pts(2) n_pts(end-1)], [0.58 1.05], '-b','linewidth',3);
axis([200 550 0.5 1.1]);
xlabel('number of timeticks (N)');
ylabel('Wall clock time (s)');
set(gca, 'fontsize', 30);

saveas(gcf, [data_out_folder dataname '_scalability'],'png');
