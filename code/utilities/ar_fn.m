function [pred_out] = ar_fn(ar_window,x, x_ori)

C = zeros(ar_window+1,size(x,2));
ar_mat = zeros(size(x,1) - ar_window, ar_window);
pred = zeros(size(x,1) - ar_window,1);
for i = 1: size(x,1) - ar_window
    ar_mat(i,:) = (x(i:i+ar_window-1,1))';
    pred(i,1) = x(i+ar_window,1);
end
ar_mat = [ones(size(ar_mat,1),1) ar_mat];

% pred = ar_mat * C
C(:,1) = ar_mat \ pred;

nn = length(x);
% now predict
for i = 1: length(x_ori) - nn
    pred_out(i,1) = [1; x_ori(nn-ar_window+i:nn+i-1)]' * C(:,1);
end