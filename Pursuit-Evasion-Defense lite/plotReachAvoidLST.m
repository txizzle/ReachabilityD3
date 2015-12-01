% [ reach, g, target, tau, obstacle] = reachAvoidSquareLST;

tPlot = -[0.05 0.15 0.2 0.4 0.5];

figure;
for i  = 1:length(tPlot)
    subplot(3,2,i)
    
    iPlot = find(tau>tPlot(i), 1 ,'first');
    [~, hr] = contour(g.xs{1}, g.xs{2}, reach(:,:,iPlot), [0 0], 'b'); hold on
    set(gca, 'fontsize', 16)
    [~, ho] = contour(g.xs{1}, g.xs{2}, obstacle(:,:,iPlot), [0 0], 'k', 'linewidth',2);
    [~, ht] = contour(g.xs{1}, g.xs{2}, target(:,:,iPlot), [0 0], 'color', [0 0.5 0], 'linewidth',2);
    
    title(['t=' num2str(tPlot(i)+0.5, '%4.2f')], 'fontsize', 14)
    axis square
end

legend([hr ho ht], {'Reach-Avoid Set', 'Physical Obstacle', 'Target Set'})
% return
% Subplot spacing
subP_size = 0.325;
subP_xmin = 0.05;
subP_ymin = 0.1;
subP_xgap = 0.03;
subP_ygap = 0.15;

subP_pos = [subP_xmin               subP_ymin+subP_size+subP_ygap       subP_size subP_size;
    subP_xmin+subP_size+subP_xgap   subP_ymin+subP_size+subP_ygap   subP_size subP_size;
    subP_xmin                       subP_ymin                      subP_size subP_size;
    subP_xmin+subP_size+subP_xgap   subP_ymin                      subP_size subP_size];

% for j = 1:4
%     subplot(3,2,j)
%     set(gca,'position',subP_pos(j,:))
% end

pos = get(gcf,'position');
set(gcf,'position',[pos(1) pos(2) 600 800]);
set(legend,'units','pixels','position',[330 100 225 150],'fontsize',14)
% legend('boxoff')