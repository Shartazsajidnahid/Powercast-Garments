function [out] = BIC(error, k, N, bic_param)

% BIC = \chai^2 + k*log2(n)

% out =  n * log( 1/n * sum((error.^2)) ) + k * log(n);

% tol = 10^(2);
% if tol < sum((error.^2))
%     error_sum = sum((error.^2));
% else
%     error_sum = tol;
% end
% 
% out =  n * log( 1/n * error_sum ) + k * log(n);



%% modified as of Sep 2016
% out = sum((error.^2)) /2/(0.0001)^2+ k * log(N); % for after - 0829 , CMU_data 0914_2
% normalization
%  out = sum((error.^2)) / 2/(0.01)^2 + k * log(N); % for time - 0906
%  script
% out = sum((error.^2)) / 2/(0.001)^2 + k * log(N); % for the latest  - 0914 script


% out = sum((error.^2))  / 2/(0.0001)^2  + k * log(N); % for after - 0829 , CMU_data 0914_2

% out = sum((error.^2))  / 2/(0.001)^2  + k * log(N); % for CMU_data 0914_2 case 1 2 flex

% out = sum((error.^2)) /2/(0.001)^2 + k * log(N); % for nmf c1

out = sum((error.^2)) / 2 / bic_param^2 + k * log(N);