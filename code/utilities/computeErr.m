function err_check_re = computeErr(I_real1, I_real_re_check)

err_check_re = sqrt(sum((I_real1-I_real_re_check).^2) / sum((I_real1).^2) );