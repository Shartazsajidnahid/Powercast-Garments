function plot_forecast_res(M_new,n_d,n_d_pred,I_real2,I_imag2,I_real_re,I_imag_re, label_tmp2,t_tot,forecast_fn, fontsize_in)

t_tot2 = length(I_real2);

% plot M new U
k = 1;
while k <= size(M_new.U{1},2)
    if rem(k,2) == 1
        figure;
    end
    subplot(2,1,rem(k-1,2)+1);hold on;
    plot(M_new.U{1}(:,k));
    plot([n_d+1:n_d+n_d_pred],M_new.U{1}(n_d+1:n_d+n_d_pred,k), 'rx');
    title(['Uhat_' num2str(k)]);
    set(gca,'XTick',[1:n_d+n_d_pred]);
    set(gca,'XTickLabel',label_tmp2);
    set(gca,'XTickLabelRotation',90);
    for i = 1: ceil(n_d/7)
        plot([4+(i-1)*7 4+(i-1)*7], [min(M_new.U{1}(:,k)) max(M_new.U{1}(:,k))], ':m');
    end
    k = k + 1;
end


figure;
subplot(2,1,1);
plot(I_real2,'k', 'linewidth',2); hold on;
plot((t_tot+1:t_tot2),I_real_re(t_tot+1:t_tot2),':r', 'linewidth',2);
set(gca,'XTick',[1:24:t_tot2]);
set(gca,'XTickLabel',label_tmp2);
set(gca,'XTickLabelRotation',90);
for i = 1: ceil(t_tot2/24)
    plot([(i-1)*24 (i-1)*24], [min(I_real2) max(I_real2)], ':k');
end
ylabel('I_r [A]', 'fontsize', fontsize_in);
xlabel('days', 'fontsize', fontsize_in);

subplot(2,1,2);
plot(I_imag2,'k', 'linewidth',2); hold on;
plot((t_tot+1:t_tot2),I_imag_re(t_tot+1:t_tot2),':r', 'linewidth',2);
set(gca,'XTick',[1:24:t_tot2]);
set(gca,'XTickLabel',label_tmp2);
set(gca,'XTickLabelRotation',90);
for i = 1: ceil(t_tot2/24)
    plot([(i-1)*24 (i-1)*24], [min(I_imag2) max(I_imag2)], ':k');
end
ylabel('I_i [A]', 'fontsize', fontsize_in);
xlabel('days', 'fontsize', fontsize_in);
legend({'Truth',forecast_fn},'orientation','horizontal','location','none', 'position', [0.10826 0.46839 0.81217 0.06990], 'fontsize', fontsize_in);
    

