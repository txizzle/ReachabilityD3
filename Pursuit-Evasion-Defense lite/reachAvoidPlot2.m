clear all; close all;
load('3+1D_tv.mat')
load('3D_tv.mat')

captureRadius = 0.1;
tMax = 1;
xd = 0.05; 
yd = -0.7;

vtarget = 1.5;
vobs = 0.5;
N2D = 200;    %   ax ay dx dy
g2D = proj2D(g, [ 0  0  1  1 ], N2D);

collision2D = cell(4,1);
initTarget.xmin = 0.4;
initTarget.xmax = 0.6;
initTarget.ymin = -0.8;
initTarget.ymax = -0.6;

target = shapeRectangleByCorners(g2D,[initTarget.xmin initTarget.ymin], [initTarget.xmax initTarget.ymax]);

initObs.xmin = -0.2;
initObs.xmax = 0;
initObs.ymin = -0.5;
initObs.ymax = 0.5;

figure;
plotnum = 0;
% for i = length(tau3D):-1:1
for i = [160 110 50 1]
    tNow = -tau3D(i);
    plotnum = plotnum+1;
    subplot(2,2,plotnum)
   
    [~, data2D] = proj2D(g3D, [0 0 1], N2D, reach3D(:,:,:,i), yd);
    [~, h3D] = contour(g2D.xs{1},g2D.xs{2}, data2D, [0 0], 'r'); hold on
    set(gca, 'fontsize', 16)
    
    plot(xd, yd, 'b*')
    collision2D{i} = sqrt((g2D.xs{1} - xd).^2 + (g2D.xs{2} - yd).^2) - captureRadius;
    [~, hc] = contour(g2D.xs{1},g2D.xs{2}, collision2D{i}, [0 0], 'k--');

    target = shapeRectangleByCorners(g2D,[initTarget.xmin initTarget.ymin+vtarget*(tMax-tNow)], [initTarget.xmax initTarget.ymax+vtarget*(tMax-tNow)]);
    [~, ht] = contour(g2D.xs{1},g2D.xs{2},target, [0 0], 'color', [0 0.5 0], 'linewidth', 2);
    
    obstacle = shapeRectangleByCorners(g2D,[initObs.xmin initObs.ymin-vobs*(tMax-tNow)], [initObs.xmax initObs.ymax]);
    [~, ho] = contour(g2D.xs{1},g2D.xs{2},obstacle,[0 0], 'k', 'linewidth', 2);
        
    title(['t=' num2str(tau3D(i)+tMax,'%4.2f')],'fontsize',14)
%     title(['i=' num2str(i)],'fontsize',14)
    axis(g.axis);
    axis square    
    drawnow
    
%     export_fig(['3Dsim2/' num2str(plotnum)], '-png')
    hold off
%     if i>1
%         delete(h3D)
%         delete(ht)
%         delete(ho)
%     end
end

% return
legend([h3D hc ht ho], {'Reach-avoid set', 'Capture set','Target set','Physical obstacle'})

% Subplot spacing
subP_size = 0.375;
subP_xmin = 0.1;
subP_ymin = 0.2;
subP_xgap = 0.1;
subP_ygap = 0;

subP_pos = [subP_xmin               subP_ymin+subP_size+subP_ygap       subP_size subP_size;
    subP_xmin+subP_size+subP_xgap   subP_ymin+subP_size+subP_ygap   subP_size subP_size;
    subP_xmin                       subP_ymin                      subP_size subP_size;
    subP_xmin+subP_size+subP_xgap   subP_ymin                      subP_size subP_size];

for j = 1:4 % Fix defender position
    subplot(2,2,j)
    set(gca,'position',subP_pos(j,:))
end

pos = get(gcf,'position');
set(gcf,'position',[pos(1) pos(2) 600 800]);
set(legend,'units','pixels','position',[150 10 225 150])
% legend('box')
% legend('boxoff')