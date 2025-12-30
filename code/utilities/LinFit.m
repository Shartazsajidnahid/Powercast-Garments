function [b, yhat, error] = LinFit(xtick, pts)


xtick2 = [ones(size(xtick,1),1) xtick];



b = xtick2 \ pts;
yhat = xtick2 * b; % linear prediction


error = pts - yhat;
