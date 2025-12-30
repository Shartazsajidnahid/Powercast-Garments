function plot_parafac_res(M,n_d, label_tmp,fontsize_in,linewidth_in, plot_save_flag,data_out_folder,dataname)


% plot lambda
figure;
bar(M.lambda);
title('lambdas');
set(gca, 'fontsize', fontsize_in);


color_all = ['k','g','b','r','c','m','y'];


% u vectors in each subplot
k = 1;
while k <= size(M.U{1},2)
    if rem(k,2) == 1
        figure;
    end
    subplot(2,1,rem(k-1,2)+1);plot(M.U{1}(:,k)); title(['U_' num2str(k)]);
    hold on;
    set(gca,'XTick',[1:n_d]);
    set(gca,'XTickLabel',label_tmp);
    set(gca,'XTickLabelRotation',90);
    set(gca, 'fontsize', fontsize_in);
    for i = 1: ceil(n_d/7)
        plot([4+(i-1)*7 4+(i-1)*7], [min(M.U{1}(:,k)) max(M.U{1}(:,k))], ':m');
    end
    k = k + 1;
end


% u vectors on top of each other
ll = cell(1);
figure; hold on;
for i = 1: size(M.U{1},2)
    plot(M.U{1}(:,i), 'linewidth',linewidth_in);
    ll{i} = ['u' num2str(i)];
end
set(gca,'XTick',[1:n_d]);
set(gca,'XTickLabel',label_tmp);
set(gca,'XTickLabelRotation',90);
set(gca, 'fontsize', fontsize_in);
legend(ll, 'fontsize', fontsize_in);
xlabel('days');
for i = 1: floor(length(M.U{1}(:,i))/7)
    plot([4+(i-1)*7 4+(i-1)*7], [min(M.U{1}(:,1)) max(M.U{1}(:,1))], ':m');
end

if plot_save_flag
    saveas(gcf, [data_out_folder dataname '_u'],'png');
end



% v vectors in each subplot
k = 1;
while k <= size(M.U{2},2)
    if rem(k,2) == 1
        figure;
    end
    subplot(2,1,rem(k-1,2)+1);plot(M.U{2}(:,k)); title(['V_' num2str(k)]);
    k = k + 1;
end



% v vectors on top of each other
ll = cell(1);
figure; hold on;
for i = 1: size(M.U{2},2)
    plot(M.U{2}(:,i), 'linewidth',linewidth_in);
    ll{i} = ['v' num2str(i)];
end
set(gca, 'fontsize', fontsize_in);
legend(ll, 'fontsize', fontsize_in);
xlabel('hour of a day');
if plot_save_flag
    saveas(gcf, [data_out_folder dataname '_v'],'png');
end



% w values for each rank component
clear xin
inin = [];
for i = 1: size(M.U{2},2)
    inin = [inin; M.U{3}(1,i) M.U{3}(3,i)];
    xin{i} = ['w' num2str(i)];
end
figure;
b = bar( inin);

set(gca, 'XTick', [1 2]);
set(gca, 'XTickLabel', xin);
set(gca, 'fontsize', fontsize_in+20);
legend({'G','B'}, 'fontsize', fontsize_in+20);

if plot_save_flag
    saveas(gcf, [data_out_folder dataname '_w'],'png');
end



BBB = M.U{3}(3,1) * M.U{2}(:,1) + M.U{3}(3,2) * M.U{2}(:,2);
GGG = M.U{3}(1,1) * M.U{2}(:,1) + M.U{3}(1,2) * M.U{2}(:,2);
figure; hold on; plot(GGG); plot(BBB);
figure;
subplot(2,1,1); plot(GGG,'k', 'linewidth',linewidth_in);
title('Dailiy G unit');
set(gca, 'fontsize', fontsize_in);
subplot(2,1,2); plot(BBB,'k',  'linewidth',linewidth_in);
title('Daily B unit');
set(gca, 'fontsize', fontsize_in);
disp('G vs B');
if plot_save_flag
    saveas(gcf, [data_out_folder dataname '_G_B'],'png');
end


figure;
hold on;
plot(GGG,'-xk', 'linewidth',linewidth_in, 'markersize',10);
plot(BBB,':k',  'linewidth',linewidth_in, 'markersize',10);
legend('Daily G unit', 'Daily B unit','location','northwest');
set(gca, 'fontsize', fontsize_in);
xlabel('Hour of the day');
if plot_save_flag
    saveas(gcf, [data_out_folder dataname '_G_B2'],'png');
end



ttt = cos(atan(BBB ./ GGG));
figure; 
plot(ttt, 'linewidth', linewidth_in);
set(gca, 'fontsize', fontsize_in);
xlabel('Hour of a day');
ylabel('Power factor');
if plot_save_flag
    saveas(gcf, [data_out_folder dataname '_G_B_ratio'],'png');
end


