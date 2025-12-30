function plot_forecast_whatif_mag_res(I_real2,I_imag2,I_real_re,I_imag_re, I_real_whatif, I_imag_whatif, label_tmp2,t_tot,forecast_fn, Gfactors, Bfactors, fontsize_in, plot_save_flag,data_out_folder, dataname,n_d_pred,tt)

t_tot2 = length(I_real2);

color_all = {'b','g','c','m'};

figure; hold on;
% h1 = plot(I_real2(1:t_tot), 'color',[0.5 0.5 0.5], 'linewidth',2); hold on;
% h2 = plot(I_real2(t_tot+1:t_tot2),'ko', 'markersize',10);
h2 = plot(I_real2(t_tot+1:t_tot2),'.k', 'markersize',20);
h3 = plot(I_real_re(t_tot+1:t_tot2),'-r', 'linewidth',2);
k = 1;
for g_idx = 1: length(Gfactors)
    for b_idx = 1: length(Bfactors)
        I_real_now = I_real_whatif{g_idx, b_idx};
        plot(I_real_now(t_tot+1:t_tot2),'color',color_all{k},'linestyle',':', 'linewidth',2);
        k = k + 1;
    end
end
axis off

if plot_save_flag
    saveas(gcf, [data_out_folder dataname '_whatif_d' num2str(n_d_pred) '_' forecast_fn '_Ir_mag' '_t' num2str(tt)],'png');
end

figure;hold on;
% subplot(2,1,2);
% h1 = plot(I_imag2(1:t_tot), 'color',[0.5 0.5 0.5], 'linewidth',2); hold on;
% h2 = plot(I_imag2(t_tot+1:t_tot2),'ko', 'markersize',10);
h2 = plot(I_imag2(t_tot+1:t_tot2),'.k', 'markersize',20);
h3 = plot(I_imag_re(t_tot+1:t_tot2),'-r', 'linewidth',2);
k = 1;
legend_entry = [h2 h3];
for g_idx = 1: length(Gfactors)
    for b_idx = 1: length(Bfactors)
        I_imag_now = I_imag_whatif{g_idx, b_idx};
        hh{k} = plot(I_imag_now(t_tot+1:t_tot2),'color',color_all{k},'linestyle',':', 'linewidth',2);
        legend_entry = [legend_entry hh{k}];
        k = k + 1;
    end
end
axis off

if plot_save_flag
    saveas(gcf, [data_out_folder dataname '_whatif_d' num2str(n_d_pred) '_' forecast_fn '_Ii_mag' '_t' num2str(tt)],'png');
end