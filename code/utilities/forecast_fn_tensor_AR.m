function [M_new] = forecast_fn_tensor_AR(M, ar_window, n_d_pred, n_d)

r = size(M.U{1},2);

C = zeros(ar_window+1,size(M.U{1},2));
for r_idx = 1: size(M.U{1},2)
    ar_mat = zeros(size(M.U{1},1) - ar_window, ar_window);
    pred = zeros(size(M.U{1},1) - ar_window,1);
    for i = 1: size(M.U{1},1) - ar_window
        ar_mat(i,:) = (M.U{1}(i:i+ar_window-1,r_idx))';
        pred(i,1) = M.U{1}(i+ar_window,r_idx);
    end
    ar_mat = [ones(size(ar_mat,1),1) ar_mat];
    
    % pred = ar_mat * C
    C(:,r_idx) = ar_mat \ pred;
end

% ar_window = 1;
% C = zeros(ar_window+1,size(M.U{1},2));
% for r_idx = 1: size(M.U{1},2)
%     ar_mat = zeros(size(M.U{1},1) - ar_window, ar_window);
%     pred = zeros(size(M.U{1},1) - ar_window,r);
%     for i = 1: size(M.U{1},1) - ar_window
%         ar_mat(i,:) = (M.U{1}(i:i+ar_window-1,:))';
%         pred(i,1:r) = M.U{1}(i+ar_window,:);
%     end
%     ar_mat = [ones(size(ar_mat,1),1) ar_mat];
%
%     % pred = ar_mat * C
%     C = ar_mat \ pred;
% end


% now predict
for r_idx = 1: r
    for i = 1: n_d_pred
        %         M_U_pred(i,r_idx) = [1; M_gt.U{1}(n_d-ar_window+i:n_d+i-1,r)]' * C(:,r_idx);
        if i == 1
            data_in1 = M.U{1}(n_d-ar_window+i:n_d+i-1,r_idx);
            data_in = data_in1;
        else
            %             data_in = [data_in(2:end); M_U_pred(i-1-ar_window+1:i-1,r_idx)];
            data_in = [data_in1; M_U_pred(1:i-1,r_idx)];
            data_in = data_in(end - ar_window+1:end);
        end
        M_U_pred(i,r_idx) = [1; data_in]' * C(:,r_idx);
    end
end


M_new.U{1} = M.U{1};
M_new.U{2} = M.U{2};
M_new.U{3} = M.U{3};
M_new.lambda = M.lambda;
M_new.U{1} = [M_new.U{1}; M_U_pred];
