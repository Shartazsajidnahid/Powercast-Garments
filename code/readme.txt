Author: Hyun Ah Song (hyunahs@cs.cmu.edu)
Date: April 21, 2017

This is source code for PowerCast algorithm in paper "PowerCast: Mining and Forecasting Power Grid Sequences", submitted to PKDD'17.

It requires:
1) Matlab Tensor Toolbox:
Brett W. Bader, Tamara G. Kolda and others. MATLAB Tensor Toolbox Version 2.6, Available online, February 2015. URL: http://www.sandia.gov/~tgkolda/TensorToolbox/. 

2) R
https://www.r-project.org/



To generate the figures in the paper:
1) download Matlab Tensor Toolbox from above url, and place it in the folder: /PowerCast/tensor_toolbox/
2) install R
3) open "/PowerCast/code/experiments/arima_seasonal_u2.R" script and change the 2nd line to be your directory where PowerCast code is
    i.e. setwd('(your directory)/PowerCast/code/tmp/')
4) open "forecast_arima_seasonal_u2.m" script and change the 7nd line to be directory where you can execute R program
    e.g. system('/usr/local/bin/R CMD BATCH arima_seasonal_u2.R');
5) to generate figures for:
    - Q1 interpretability and Q2 forecasting: run "/PowerCast/code/experiments/PowerCast.m"
    - Q2 forecasting under what-if scenarios: run "/PowerCast/code/experiments/whatif.m"
    - Q3 anomaly detection: run "/PowerCast/code/experiments/anomaly_detection.m"
    - Q4 anomaly explanation: run "/PowerCast/code/experiments/anomaly_whatif.m"
    - Scalability: run "/PowerCast/code/experiments/scalability.m"

All figures will be saved in "./PowerCast/plots" folder.

