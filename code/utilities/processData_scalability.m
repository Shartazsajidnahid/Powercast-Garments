function [V_real1, I_real1, V_imag1, I_imag1, myTensor, V_real2, I_real2, V_imag2, I_imag2, n_d, label_tmp_unit, label_tmp, label_tmp2, t_d] = processData_scalability(...
    dataname, n_d_pred, n_d_removed, sigma, window_size, fontsize_in, plot_save_flag, linewidth_in,data_out_folder)


% set the starting & ending point for each dataset
if strcmp(dataname, 'CMU')
    
    Dir = '../data/CMU_data/';
    Datacode = 'int3';
    Dataname = sprintf('%s_data.mat', Datacode); % I_real, I_imag, V_real, V_imag
    % Dataname = 'data_02222017_int3_data.mat';
    load([Dir Dataname]);
    Data_time_name = sprintf('%s_dates.mat', Datacode);
    load([Dir Data_time_name]); % timestamp
    
    
    start_pt = 13%10%13;
    %     n_d_removed = 0
    end_pt = 660 - 24*n_d_removed - 24*n_d_pred%%%%660 - 24*n_d_pred;%657 - 24*n_d_pred%660 - 24*n_d_pred;% 660;
    
    
    % input data
    V_real1 = V_real(start_pt:end_pt);
    I_real1 = I_real(start_pt:end_pt);
    V_imag1 = V_imag(start_pt:end_pt);
    I_imag1 = I_imag(start_pt:end_pt);
    
    
    
    t_tot = length(V_real1);
    t_d = 24; % hourly basis: 24 hours/day
    n_d = t_tot/t_d; % hourly basis

    
    
    %%%%% ground-truth
    start_pt = 13%10%13;
    end_pt = 660 - 24*n_d_removed%%%%660%657%660;% 660; % V_real is weird - too low for the last two days! skip these!
    
    V_real2 = V_real(start_pt:end_pt);
    I_real2 = I_real(start_pt:end_pt);
    V_imag2 = V_imag(start_pt:end_pt);
    I_imag2 = I_imag(start_pt:end_pt);
    %%%%%%
    
    
    
    
    
elseif strcmp(dataname, 'LBNL')
    
    
    fnum = 120;
    Dir = '../data/LBNL/';
    Datacode = 'lbnl_bank';
    Dataname = sprintf('bank%d.mat', fnum);
    seg_method = 'greedy';
    downsample_factor = 3600;
    points_per_day = 120 * 60 * 60 * 24 / fnum / downsample_factor; % 120 Hz
    
    
    % BIC_params = 2.^(-2:.5:3);
    BIC_params = 2^3;
    
    num_trials = length(BIC_params);
    
    load([Dir Dataname]);
    Vi = V_imag; Vr = V_real; Ii = I_imag; Ir = I_real;
    
    seqlen = length(Vr) - mod(length(Vr), downsample_factor);
    V_real = mean(reshape(Vr(1:seqlen), downsample_factor, []), 1)';
    V_imag = mean(reshape(Vi(1:seqlen), downsample_factor, []), 1)';
    I_real = mean(reshape(Ir(1:seqlen), downsample_factor, []), 1)';
    I_imag = mean(reshape(Ii(1:seqlen), downsample_factor, []), 1)';
    %%
    
    start_pt = 1;
    n_d_ori = floor(length(V_real)/24);
    end_pt = n_d_ori*24 - 24*n_d_removed - 24*n_d_pred%%%%660 - 24*n_d_pred;%657 - 24*n_d_pred%660 - 24*n_d_pred;% 660;
    
    
    % input data
    V_real1 = V_real(start_pt:end_pt);
    I_real1 = I_real(start_pt:end_pt);
    V_imag1 = V_imag(start_pt:end_pt);
    I_imag1 = I_imag(start_pt:end_pt);
    
    
    
    t_tot = length(V_real1);
    t_d = 24; % hourly basis: 24 hours/day
    n_d = t_tot/t_d; % hourly basis
    
    
    % Thu, 01 Oct 2015 00:00:00 GMT, and the timestamps are in 1 second intervals
    
    % timestamp for data1
    label_tmp_unit = {'T' 'F' 'S' 'S' 'M' 'T' 'W'};
    if n_d / 7 > 1
        label_tmp = repmat(label_tmp_unit, 1,floor(n_d/7));
        label_tmp = [label_tmp label_tmp{1:n_d-size(label_tmp,2)}];
    else
        label_tmp = label_tmp_unit(1:n_d);
    end
    
    unix_time = 1443657600;
    for i = 1:n_d
        datenow = datestr((unix_time + 24*(i-1)*60*60)/86400 + datenum(1970,1,1));
        timestamp_days_sel(i,1) =  strcat(datenow,'-',label_tmp(i));
    end
    label_tmp = timestamp_days_sel;
    
    % timestamp for data2
    n_d2 = n_d + n_d_pred;
    if n_d2 / 7 > 1
        label_tmp2 = repmat(label_tmp_unit, 1,floor(n_d2/7));
        label_tmp2 = [label_tmp2 label_tmp2{1:n_d2-size(label_tmp2,2)}];
    else
        label_tmp2 = label_tmp_unit(1:n_d2);
    end
    
    for i = 1:n_d+n_d_pred
        datenow = datestr((unix_time + 24*(i-1)*60*60)/86400 + datenum(1970,1,1));
        timestamp_days_sel(i,1) =  strcat(datenow,'-',label_tmp2(i));
    end
    label_tmp2 = timestamp_days_sel;
    
    
    
    
    
    
    %%%%% ground-truth
    start_pt = 1%10%13;
    end_pt = n_d_ori*24 - 24*n_d_removed%%%%660%657%660;% 660; % V_real is weird - too low for the last two days! skip these!
    
    V_real2 = V_real(start_pt:end_pt);
    I_real2 = I_real(start_pt:end_pt);
    V_imag2 = V_imag(start_pt:end_pt);
    I_imag2 = I_imag(start_pt:end_pt);
    %%%%%%
    
    
    
    
    
end


data_proc_method = 1;
% window_size = 5; % should be odd number
% sigma = 0.5;

myTensor = construct_tensor(V_real1,I_real1,V_imag1,I_imag1,n_d,t_d,data_proc_method,window_size,sigma, dataname);



% let's plot the parameters - tensor values
G = myTensor(:,:,1);
alpha_r = myTensor(:,:,2);
B = myTensor(:,:,3);
alpha_i = myTensor(:,:,4);

G_mat = reshape(G',1,size(G,1)*size(G,2));
alpha_r_mat = reshape(alpha_r',1,size(alpha_r,1)*size(alpha_r,2));
B_mat = reshape(B',1,size(B,1)*size(B,2));
alpha_i_mat = reshape(alpha_i',1,size(alpha_i,1)*size(alpha_i,2));

