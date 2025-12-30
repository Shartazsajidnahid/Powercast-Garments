function [V_real1, I_real1, V_imag1, I_imag1, myTensor, V_real2, I_real2, V_imag2, I_imag2, n_d, label_tmp_unit, label_tmp, label_tmp2, t_d] = processData(...
    dataname, n_d_pred, n_d_removed, sigma, window_size, fontsize_in, plot_save_flag, linewidth_in,data_out_folder)


% set the starting & ending point for each dataset
if strcmp(dataname, 'CMU')
    
    Dir = '../../data/CMU_data/';
    Datacode = 'int3';
    Dataname = sprintf('%s_data.mat', Datacode); % I_real, I_imag, V_real, V_imag
    % Dataname = 'data_02222017_int3_data.mat';
    load([Dir Dataname]);
    Data_time_name = sprintf('%s_dates.mat', Datacode);
    load([Dir Data_time_name]); % timestamp
    for i = 1: size(timestamp,1)
        timestamp_fin{i,1} = timestamp{i,1}(1:10);
    end
    timestamp_days = unique(timestamp_fin);
    
    
    start_pt = 13;%10%13;
    %     n_d_removed = 0
    end_pt = 660 - 24*n_d_removed - 24*n_d_pred;%%%%660 - 24*n_d_pred;%657 - 24*n_d_pred%660 - 24*n_d_pred;% 660;
    
    
    % input data
    V_real1 = V_real(start_pt:end_pt);
    I_real1 = I_real(start_pt:end_pt);
    V_imag1 = V_imag(start_pt:end_pt);
    I_imag1 = I_imag(start_pt:end_pt);
    
    
    
    t_tot = length(V_real1);
    t_d = 24; % hourly basis: 24 hours/day
    n_d = t_tot/t_d; % hourly basis
    
    
    
    
    
    % timestamp for data1
    label_tmp_unit = {'F' 'S' 'S' 'M' 'T' 'W' 'T'};
    label_tmp = repmat(label_tmp_unit, 1,floor(n_d/7));
    label_tmp = [label_tmp label_tmp{1:n_d-size(label_tmp,2)}];
    
    k = 1;
    for i = 2: size(timestamp_days,1)-1 -n_d_removed - n_d_pred
        timestamp_days_sel(k,1) =  strcat(timestamp_days{i,1},'-',label_tmp(k));
        k = k + 1;
    end
    label_tmp = timestamp_days_sel;
    
    %     label_tmp_unit2 = {'Fri' 'Sat' 'Sun' 'Mon' 'Tue' 'Wed' 'Thu'};
    %     label_tmp2 = repmat(label_tmp_unit2, 1,floor(n_d/7));
    %     label_tmp2 = [label_tmp2 label_tmp2{1:n_d-size(label_tmp2,2)}];
    
    % timestamp for data2
    n_d2 = n_d + n_d_pred;
    label_tmp2 = repmat(label_tmp_unit, 1,floor(n_d2/7));
    label_tmp2 = [label_tmp2 label_tmp2{1:n_d2-size(label_tmp2,2)}];
    
    k = 1;
    for i = 2: size(timestamp_days,1)-1 - n_d_removed
        timestamp_days_sel(k,1) =  strcat(timestamp_days{i,1},'-',label_tmp2(k));
        k = k + 1;
    end
    label_tmp2 = timestamp_days_sel;
    
    
    
    
    %%%%% ground-truth
    start_pt = 13;%10%13;
    end_pt = 660 - 24*n_d_removed;%%%%660%657%660;% 660; % V_real is weird - too low for the last two days! skip these!
    
    V_real2 = V_real(start_pt:end_pt);
    I_real2 = I_real(start_pt:end_pt);
    V_imag2 = V_imag(start_pt:end_pt);
    I_imag2 = I_imag(start_pt:end_pt);
    %%%%%%
    

elseif strcmp(dataname, 'Garments')
    
    Dir = '../../data/Garments/';
    Datacode = 'Garments';
    Dataname = sprintf('%s.mat', Datacode); % I_real, I_imag, V_real, V_imag
    % Dataname = 'data_02222017_int3_data.mat';
    garments_data = load([Dir Dataname]);
    garments_data = garments_data.data;
    garments_data = garments_data(1:336, :);


    Data_time_name = 'int3_dates.mat';
    load([Dir Data_time_name]); % timestamp
    for i = 1: size(timestamp,1)
        timestamp_fin{i,1} = timestamp{i,1}(1:10);
    end
    timestamp_days = unique(timestamp_fin);

    timestamp_days = timestamp_days(1:14);
    
    
    start_pt = 1;%10%13;
    n_d_removed = 3;
    % end_pt = 330 - 24*n_d_removed - 24*n_d_pred;%%%%660 - 24*n_d_pred;%657 - 24*n_d_pred%660 - 24*n_d_pred;% 660;
    end_pt = 250;%%%%660 - 24*n_d_pred;%657 - 24*n_d_pred%660 - 24*n_d_pred;% 660;
    
    
    % garments columns
    V_real = garments_data(:,1);
    I_real = garments_data(:,2);
    V_imag = garments_data(:,3);
    I_imag = garments_data(:,4);


    % input data
    V_real1 = V_real(start_pt:end_pt);
    I_real1 = I_real(start_pt:end_pt);
    V_imag1 = V_imag(start_pt:end_pt);
    I_imag1 = I_imag(start_pt:end_pt);
    
    
    
    t_tot = length(V_real1);
    t_d = 10; % hourly basis: 24 hours/day
    n_d = t_tot/t_d; % hourly basis
    
    
    
    
    
    % timestamp for data1
    label_tmp_unit = {'F' 'S' 'S' 'M' 'T' 'W' 'T'};
    label_tmp = repmat(label_tmp_unit, 1,floor(n_d/7));
    label_tmp = [label_tmp label_tmp{1:n_d-size(label_tmp,2)}];
    
    k = 1;
    for i = 2: size(timestamp_days,1)-1 -n_d_removed - n_d_pred
        timestamp_days_sel(k,1) =  strcat(timestamp_days{i,1},'-',label_tmp(k));
        k = k + 1;
    end
    label_tmp = timestamp_days_sel;
    
    %     label_tmp_unit2 = {'Fri' 'Sat' 'Sun' 'Mon' 'Tue' 'Wed' 'Thu'};
    %     label_tmp2 = repmat(label_tmp_unit2, 1,floor(n_d/7));
    %     label_tmp2 = [label_tmp2 label_tmp2{1:n_d-size(label_tmp2,2)}];
    
    % timestamp for data2
    n_d2 = n_d + n_d_pred;
    label_tmp2 = repmat(label_tmp_unit, 1,floor(n_d2/7));
    label_tmp2 = [label_tmp2 label_tmp2{1:n_d2-size(label_tmp2,2)}];
    
    k = 1;
    for i = 2: size(timestamp_days,1)-1 - n_d_removed
        timestamp_days_sel(k,1) =  strcat(timestamp_days{i,1},'-',label_tmp2(k));
        k = k + 1;
    end
    label_tmp2 = timestamp_days_sel;
    
    
    
    
    %%%%% ground-truth
    start_pt = 1;%10%13;
    % end_pt = 660 - 24*n_d_removed;%%%%660%657%660;% 660; % V_real is weird - too low for the last two days! skip these!
    end_pt = 310;
    
    V_real2 = V_real(start_pt:end_pt);
    I_real2 = I_real(start_pt:end_pt);
    V_imag2 = V_imag(start_pt:end_pt);
    I_imag2 = I_imag(start_pt:end_pt);
    %%%%%%
    
    
    
elseif strcmp(dataname, 'LBNL')
    
    
    fnum = 120;
    Dir = '../../data/LBNL/';
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
    end_pt = n_d_ori*24 - 24*n_d_removed - 24*n_d_pred;%%%%660 - 24*n_d_pred;%657 - 24*n_d_pred%660 - 24*n_d_pred;% 660;;
    
    
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
    start_pt = 1;%10%13;
    end_pt = n_d_ori*24 - 24*n_d_removed;%%%%660%657%660;% 660; % V_real is weird - too low for the last two days! skip these!
    
    V_real2 = V_real(start_pt:end_pt);
    I_real2 = I_real(start_pt:end_pt);
    V_imag2 = V_imag(start_pt:end_pt);
    I_imag2 = I_imag(start_pt:end_pt);
    %%%%%%
    
    
    
    
    
end


data_proc_method = 1;
window_size = 3; % should be odd number
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

figure;
subplot(4,1,1); plot(G_mat, 'linewidth',linewidth_in);hold on; title('G');
for i = 1: ceil(length(G_mat)/24)
    plot([0+(i-1)*24 0+(i-1)*24], [min(G_mat) max(G_mat)], ':m');
end
set(gca, 'fontsize', fontsize_in);
subplot(4,1,2); plot(alpha_r_mat, 'linewidth',linewidth_in);hold on;title('alpha_r');
for i = 1: ceil(length(G_mat)/24)
    plot([0+(i-1)*24 0+(i-1)*24], [min(alpha_r_mat) max(alpha_r_mat)], ':m');
end
set(gca, 'fontsize', fontsize_in);
subplot(4,1,3); plot(B_mat, 'linewidth',linewidth_in);hold on;title('B');
for i = 1: ceil(length(G_mat)/24)
    plot([0+(i-1)*24 0+(i-1)*24], [min(B_mat) max(B_mat)], ':m');
end
set(gca, 'fontsize', fontsize_in);
subplot(4,1,4); plot(alpha_i_mat, 'linewidth',linewidth_in); hold on;title('alpha_i');
for i = 1: ceil(length(G_mat)/24)
    plot([0+(i-1)*24 0+(i-1)*24], [min(alpha_i_mat) max(alpha_i_mat)], ':m');
end
set(gca, 'fontsize', fontsize_in);


if plot_save_flag
    saveas(gcf, [data_out_folder dataname '_parameters'],'png');
end