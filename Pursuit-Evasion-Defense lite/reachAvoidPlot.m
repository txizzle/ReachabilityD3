clear all; close all
load('3+1D_tv.mat')
load('3D_tv.mat')

captureRadius = 0.1;

xd = 0.05; 
yd = [-0.7 -0.3 0.1 0.5];

N2D = 200;
g2D = proj2D(g, [0 0 1 1], N2D);

collision2D = cell(4,1);
initTarget.xmin = 0.4;
initTarget.xmax = 0.6;
initTarget.ymin = -0.8;
initTarget.ymax = -0.6;

target = shapeRectangleByCorners(g2D,[initTarget.xmin initTarget.ymin], [initTarget.xmax initTarget.ymax]);
vtarget = 1.5;    % velocity of the target set (moving downwards in forward time)


initObs.xmin = -0.2;
initObs.xmax = 0;
initObs.ymin = -0.5;
initObs.ymax = 0.5;
obstacle = shapeRectangleByCorners(g2D,[initObs.xmin initObs.ymin], [initObs.xmax initObs.ymax]);
vobs = 0.5;       % velocity of obstacle (moving left in forward time)

tMax = 1;                  % End time.
tNow = -tau3D(1);
figure;
for i = 1:length(yd)
    subplot(2,2,i)
    [~, data2D] = proj2D(g, [0 0 1 1], N2D, reach, [yd(i) 0]);
    [~, h4D{i}] = contour(g2D.xs{1},g2D.xs{2}, data2D, [0 0], 'b:','linewidth',2);  hold on
    set(gca, 'fontsize', 16)
    
    [~, data2D] = proj2D(g3D, [0 0 1], N2D, reach3D(:,:,:,1), yd(i));
    [~, h3D{i}] = contour(g2D.xs{1},g2D.xs{2}, data2D, [0 0], 'r');

    plot(xd, yd(i), 'b*')
    collision2D{i} = sqrt((g2D.xs{1} - xd).^2 + (g2D.xs{2} - yd(i)).^2) - captureRadius;
    [~, hc{i}] = contour(g2D.xs{1},g2D.xs{2}, collision2D{i}, [0 0], 'k--');

    target = shapeRectangleByCorners(g2D,[initTarget.xmin initTarget.ymin+vtarget*(tMax-tNow)], [initTarget.xmax initTarget.ymax+vtarget*(tMax-tNow)]);
    [~, ht{i}] = contour(g2D.xs{1},g2D.xs{2},target, [0 0], 'color', [0 0.5 0], 'linewidth', 2);
    
    obstacle = shapeRectangleByCorners(g2D, [initObs.xmin initObs.ymin-vobs*(tMax-tNow)], [initObs.xmax initObs.ymax]);    
    [~, ho{i}] = contour(g2D.xs{1},g2D.xs{2},obstacle,[0 0], 'k', 'linewidth', 2);
        
    title('t=0','fontsize',14)
    axis(g.axis);
    axis square    
end

legend([h4D{4} h3D{4} hc{4} ht{4} ho{4}], {'Reach-avoid set (4D)','Reach-avoid set (3D)', 'Capture set','Target set','Physical obstacle'})

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
set(gcf,'position',[200 200 600 800]);
set(legend,'units','pixels','position',[150 10 225 150])
% legend('boxoff')