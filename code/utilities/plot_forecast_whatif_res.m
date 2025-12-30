function plot_forecast_whatif_res(I_real2,I_imag2,I_real_re,I_imag_re, I_real_whatif, I_imag_whatif, label_tmp2,t_tot,forecast_fn, Gfactors, Bfactors, fontsize_in)

t_tot2 = length(I_real2);

color_all = {'b','g','c','m'};

figure;
subplot(2,1,1);
h1 = plot(I_real2(1:t_tot), 'color',[0.5 0.5 0.5], 'linewidth',2); hold on;
% h2 = plot((t_tot+1:t_tot2),I_real2(t_tot+1:t_tot2),'o','color',[0.5 0.5 0.5], 'markersize',3);
h2 = plot((t_tot+1:t_tot2),I_real2(t_tot+1:t_tot2),'.k', 'markersize',5);
h3 = plot((t_tot+1:t_tot2),I_real_re(t_tot+1:t_tot2),'-r', 'linewidth',2);
k = 1;
for g_idx = 1: length(Gfactors)
    for b_idx = 1: length(Bfactors)
        I_real_now = I_real_whatif{g_idx, b_idx};
        plot((t_tot+1:t_tot2),I_real_now(t_tot+1:t_tot2),'color',color_all{k},'linestyle',':', 'linewidth',2);
        k = k + 1;
    end
end
set(gca,'XTick',[1:24:t_tot2]);
set(gca,'XTickLabel',label_tmp2);
set(gca,'XTickLabelRotation',90);
for i = 1: ceil(t_tot2/24)
    plot([(i-1)*24 (i-1)*24], [min(I_real2) max(I_real2)], ':k');
end
ylabel('I_r [A]', 'fontsize', fontsize_in);
xlabel('days', 'fontsize', fontsize_in);



subplot(2,1,2);
h1 = plot(I_imag2(1:t_tot), 'color',[0.5 0.5 0.5], 'linewidth',2); hold on;
% h2 = plot((t_tot+1:t_tot2),I_imag2(t_tot+1:t_tot2),'o','color',[0.5 0.5 0.5], 'markersize',3);
h2 = plot((t_tot+1:t_tot2),I_imag2(t_tot+1:t_tot2),'.k', 'markersize',5);
h3 = plot((t_tot+1:t_tot2),I_imag_re(t_tot+1:t_tot2),'-r', 'linewidth',2);
k = 1;
legend_entry = [h2 h3];
for g_idx = 1: length(Gfactors)
    for b_idx = 1: length(Bfactors)
        I_imag_now = I_imag_whatif{g_idx, b_idx};
        hh{k} = plot((t_tot+1:t_tot2),I_imag_now(t_tot+1:t_tot2),'color',color_all{k},'linestyle',':', 'linewidth',2);
        legend_entry = [legend_entry hh{k}];
        k = k + 1;
    end
end


set(gca,'XTick',[1:24:t_tot2]);
set(gca,'XTickLabel',label_tmp2);
set(gca,'XTickLabelRotation',90);
for i = 1: ceil(t_tot2/24)
    plot([(i-1)*24 (i-1)*24], [min(I_imag2) max(I_imag2)], ':k');
end
ylabel('I_i [A]', 'fontsize', fontsize_in);
xlabel('days', 'fontsize', fontsize_in);
legend_tmp = {'Truth', 'forecast'};
k = size(legend_tmp,2)+1;
for g_idx = 1: length(Gfactors)
    for b_idx = 1: length(Bfactors)
        %        legend_tmp{k} = ['G x ' num2str(Gfactors(g_idx)) ', B x ' num2str(Bfactors(b_idx))];
        if length(Gfactors)~=1
            legend_tmp{k} = ['G x ' num2str(Gfactors(g_idx))];
        elseif length(Bfactors)~=1
            legend_tmp{k} = ['B x ' num2str(Bfactors(b_idx))];
        end
        k = k + 1;
    end
end


legend(legend_entry,legend_tmp,'orientation','horizontal','location','none', 'position', [0.1622092147 0.45168 0.715 0.0462], 'fontsize', fontsize_in);

