function [error_sq, rel_error] = CompErr(xss_all, yss_all, yhat_all)

N = size(xss_all{1},1);
% compute error
for r = 1: size(xss_all,2)
    
    y = yss_all{r};
    
    yhat = [];
    for i = 1:size(yhat_all{r},2)
        yhat = [yhat; yhat_all{r}{i}];
    end
    
%     error_abs(r) = sum(abs(y-yhat))/N;
    error_sq(r) = sqrt(sum((y-yhat).^2));
    rel_error(r) = error_sq(r) / sqrt(N);
end