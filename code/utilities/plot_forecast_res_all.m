function plot_forecast_res_all(n_ds, predictions_re, predictions_imag, I_real2, I_imag2, label_tmp2, forecast_fns, dataname,data_out_folder,plot_save_flag, fontsize_in)




t_tot2 = length(I_real2);



for n_d_pred_idx = 1: length(n_ds)
    
    n_d_pred = n_ds(n_d_pred_idx);
    t_tot = t_tot2 - 24 * n_d_pred;
    
    
    
    figure;hold on;
    
    subplot(2,1,1); hold on;
    
    f1 = [1 2 3 4];
    v1 = [t_tot+1 min(I_real2); t_tot+1 max(I_real2) ;t_tot2 max(I_real2);t_tot2 min(I_real2)];
    
    patch('Faces',f1, 'Vertices',v1, 'EdgeColor','none','FaceColor',[.9 .9 .9]);
    
    title([num2str(24*n_d_pred) ' points forecasting' ' (' num2str(n_d_pred) 'days)'], 'fontsize',fontsize_in);
    h1 = plot(I_real2(1:t_tot), 'color',[0.5 0.5 0.5], 'linewidth',2);
    h2 = plot((t_tot+1:t_tot2),I_real2(t_tot+1:t_tot2),'.k','markersize',8);
    h3 = plot((t_tot+1:t_tot2),predictions_re{1,n_d_pred}(t_tot+1:t_tot2),'-r', 'linewidth',2);
    h4 = plot((t_tot+1:t_tot2),predictions_re{2,n_d_pred}(t_tot+1:t_tot2),'-m', 'linewidth',2);
    h5 = plot([t_tot+1:t_tot2], predictions_re{3,n_d_pred}, '-', 'color',[0 0 0.5],'linewidth',2);
    h6 = plot([t_tot+1:t_tot2], predictions_re{4,n_d_pred}, '-','color',[0 0.5 0], 'linewidth',2);
    
    
    set(gca,'XTick',[1:24:t_tot2]);
    set(gca,'XTickLabel',label_tmp2);
    set(gca,'XTickLabelRotation',90);
    
    for i = 1: ceil(t_tot2/24)
        plot([(i-1)*24 (i-1)*24], [min(I_real2) max(I_real2)], ':k');
    end
    ylabel('I_r [A]', 'fontsize',fontsize_in);
    xlabel('days', 'fontsize',fontsize_in);
    
    
    subplot(2,1,2);hold on;
    f1 = [1 2 3 4];
    v1 = [t_tot+1 min(I_imag2); t_tot+1 max(I_imag2) ;t_tot2 max(I_imag2);t_tot2 min(I_imag2)];
    
    patch('Faces',f1, 'Vertices',v1, 'EdgeColor','none','FaceColor',[0 0 0.1],'FaceAlpha',.05);
    
    h1 = plot(I_imag2(1:t_tot), 'color',[0.5 0.5 0.5], 'linewidth',2);
    h2 = plot((t_tot+1:t_tot2),I_imag2(t_tot+1:t_tot2),'.k', 'markersize',8);
    h3 = plot((t_tot+1:t_tot2),predictions_imag{1,n_d_pred}(t_tot+1:t_tot2),'-r', 'linewidth',2);
    h4 = plot((t_tot+1:t_tot2),predictions_imag{2,n_d_pred}(t_tot+1:t_tot2),'-m', 'linewidth',2);
    h5 = plot([t_tot+1:t_tot2], predictions_imag{3,n_d_pred}, '-', 'color',[0 0 0.5],'linewidth',2);
    h6 = plot([t_tot+1:t_tot2], predictions_imag{4,n_d_pred}, '-','color',[0 0.5 0], 'linewidth',2);
    
    
    %     plot([t_tot t_tot], [min(I_imag2) max(I_imag2)], 'k', 'linewidth',2);
    
    
    set(gca,'XTick',[1:24:t_tot2]);
    set(gca,'XTickLabel',label_tmp2);
    set(gca,'XTickLabelRotation',90);
    for i = 1: ceil(t_tot2/24)
        plot([(i-1)*24 (i-1)*24], [min(I_imag2) max(I_imag2)], ':k');
    end
    ylabel('I_i [A]', 'fontsize',fontsize_in);
    xlabel('days', 'fontsize',fontsize_in);
    
    legend([h2 h3 h4 h5 h6],{'Truth','PowerCast','PowerCast-S','AR','SAR'},'orientation','horizontal','location','none', 'position', [0.10826 0.46839 0.81217 0.06990], 'fontsize',fontsize_in);
    
    if plot_save_flag
        saveas(gcf, [data_out_folder dataname '_predictions_all_' num2str(n_d_pred) '_res'],'png');
    end
    
    
    
    
    
    % magnifying window
    figure; hold on;set(gca,'Color',[.9 .9 .9]);
    h2 = plot(I_real2(t_tot+1:t_tot2),'.k', 'linewidth',2,'markersize',40);
    h3 = plot(predictions_re{1,n_d_pred}(t_tot+1:t_tot2),'-r', 'linewidth',2);
    h4 = plot(predictions_re{2,n_d_pred}(t_tot+1:t_tot2),'-m', 'linewidth',2);
    h5 = plot(predictions_re{3,n_d_pred}, '-','color',[0 0 0.5], 'linewidth',2);
    h6 = plot(predictions_re{4,n_d_pred}, '-','color',[0 0.5 0],'linewidth',2);
    axis off
    
    if plot_save_flag
        saveas(gcf, [data_out_folder dataname '_predictions_all_magnified_re_' num2str(n_d_pred) '_res'],'png');
    end
    
    % magnifying window
    figure; hold on;set(gca,'Color',[.9 .9 .9]);
    h2 = plot(I_imag2(t_tot+1:t_tot2),'.k', 'linewidth',2,'markersize',40);
    h3 = plot(predictions_imag{1,n_d_pred}(t_tot+1:t_tot2),'-r', 'linewidth',2);
    h4 = plot(predictions_imag{2,n_d_pred}(t_tot+1:t_tot2),'-m', 'linewidth',2);
    h5 = plot(predictions_imag{3,n_d_pred}, '-','color',[0 0 0.5], 'linewidth',2);
    h6 = plot(predictions_imag{4,n_d_pred}, '-','color',[0 0.5 0],'linewidth',2);
    axis off
    
    if plot_save_flag
        saveas(gcf, [data_out_folder dataname '_predictions_all_magnified_imag_' num2str(n_d_pred) '_res'],'png');
    end
end