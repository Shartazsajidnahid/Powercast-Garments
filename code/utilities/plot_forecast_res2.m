function plot_forecast_res2(nn, pred_re, pred_imag, I_real2, I_imag2, forecast_fn,label_tmp2, fontsize_in)

t_tot2 = length(I_real2);

figure;
subplot(2,1,1); plot(I_real2, 'k', 'linewidth',2); hold on;
plot([nn+1:length(I_real2)], pred_re, ':r', 'linewidth',2);

set(gca,'XTick',[1:24:t_tot2]);
set(gca,'XTickLabel',label_tmp2);
set(gca,'XTickLabelRotation',90);
for i = 1: ceil(t_tot2/24)
    plot([(i-1)*24 (i-1)*24], [min(I_real2) max(I_real2)], ':k');
end


ylabel('I_r [A]', 'fontsize', fontsize_in);
xlabel('days', 'fontsize', fontsize_in);

subplot(2,1,2); plot(I_imag2, 'k', 'linewidth',2); hold on;
plot([nn+1:length(I_imag2)], pred_imag, ':r', 'linewidth',2);

set(gca,'XTick',[1:24:t_tot2]);
set(gca,'XTickLabel',label_tmp2);
set(gca,'XTickLabelRotation',90);
for i = 1: ceil(t_tot2/24)
    plot([(i-1)*24 (i-1)*24], [min(I_imag2) max(I_imag2)], ':k');
end

ylabel('I_i [A]', 'fontsize', fontsize_in);
xlabel('days', 'fontsize', fontsize_in);
legend({'Truth',forecast_fn},'orientation','horizontal','location','none', 'position', [0.10826 0.46839 0.81217 0.06990], 'fontsize', fontsize_in);


