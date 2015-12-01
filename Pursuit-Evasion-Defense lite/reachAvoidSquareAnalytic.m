% % function reachAvoidSquareAnalytic
clear all; close all
% Nx = [51 101 151 201 251 301]';
% reach = cell(size(Nx));
% g = cell(size(Nx));
% target = cell(size(Nx));
% tau = cell(size(Nx));
% obstacle = cell(size(Nx));
% 
% sd_reach = cell(size(Nx));
% 
% for i = 1:length(Nx)
%     [ reach{i}, g{i}, target{i}, tau{i}, obstacle{i}] = reachAvoidSquareLST(Nx(i),'veryHigh');
%     tic; sd_reach{i} = signedDistanceIterative(g{i}, reach{i}(:,:,1), 'veryHigh', 1); toc
% end
% 
% save('analytic_data.mat')
% 
% return

load('analytic_data.mat')

figure;
contour(g{end}.xs{1}, g{end}.xs{2}, sd_reach{end}, [0 0], 'b'); hold on
contour(g{end}.xs{1}, g{end}.xs{2}, target{end}(:,:,1), [0 0], 'linecolor', [0 0.5 0], 'linewidth',2)
contour(g{end}.xs{1}, g{end}.xs{2}, obstacle{end}(:,:,1), [0 0], 'k', 'linewidth',2)
axis square
xlabel('x')
ylabel('y')

velocity = 0.5;   % velocity of vehicle
vtarget = 1.5;    % velocity of the target set (moving downwards in forward time)
vobs = 1;       % velocity of obstacle (moving left in forward time)

tMax = 0.5;                  % End time.
n = 10000;       % Number of analytic points per 1 units in segment length

initTarget.xmin = -0.2;
initTarget.xmax = 0.2;
initTarget.ymin = -0.2;
initTarget.ymax = 0.2;

initObs.xmin = -0.1;
initObs.xmax = 0.1;
initObs.ymin = -0.6;
initObs.ymax = -0.4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bottom curved segment
% x = (-0.1:ds:0.1)'; % LEFT AND RIGHT SIDE
yi = initObs.ymin + vobs*tMax;
xL = initObs.xmin;
xR = initObs.xmax;

C = tMax*velocity- (initTarget.ymin-yi);

vec = [1 -1; vobs velocity]\[C-yi; velocity*yi];

d = vec(1);
yf = vec(2);

x = (-0.1:0.001:0)';
y = yf - sqrt(d^2 - min([(x-xL) (xR-x)],[],2).^2);

plength = sum(sqrt(diff(x).^2 + diff(y).^2));   % length of segment
% ds = plength/n;                              % want 1000 points per 0.1 segment length
x = linspace(-0.1,0,round(plength*n))';
y = yf - sqrt(d^2 - min([(x-xL) (xR-x)],[],2).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bottom straight segments
% x0 = [initTarget.xmin:ds:xL xR:ds:initTarget.xmax]; % LEFT AND RIGHT SIDE
plength = xL - initTarget.xmin;
% x0 = initTarget.xmin:ds:xL;
x0 = linspace(initTarget.xmin,xL,round(plength*n));
y0 = initTarget.ymin - tMax*velocity*ones(size(x0));
x = cat(1,x0',x);
y = cat(1,y0',y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Side straight segments
plength = initTarget.ymax - initTarget.ymin;
% y0 = initTarget.ymin:ds:initTarget.ymax;
y0 = linspace(initTarget.ymin,initTarget.ymax,round(plength*n));
x0 = initTarget.xmin - tMax*velocity*ones(size(y0));

x = cat(1,x0',x);
y = cat(1,y0',y);

% RIGHT SIDE
% x0 = initTarget.xmax + tMax*velocity*ones(size(y0));
% x = cat(1,x0',x);
% y = cat(1,y0',y);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lower circular segment at lower left
% Initial curve with arbitrary resolution
cc = [initTarget.xmin initTarget.ymin] ;
rr = tMax*velocity;
th = pi:0.01:3*pi/2;
xy = [ones(size(th'))*cc(1) ones(size(th'))*cc(2)] + rr * [cos(th') sin(th')];
x0 = xy(:,1);
y0 = xy(:,2);

% Refine resolution
plength = sum(sqrt(diff(x0).^2 + diff(y0).^2));   % length of segment
% ds = plength/n;                              % want 1000 points per 0.1 segment length

th = linspace(pi,3*pi/2,round(plength*n));
xy = [ones(size(th'))*cc(1) ones(size(th'))*cc(2)] + rr * [cos(th') sin(th')];
x0 = xy(:,1);
y0 = xy(:,2);

x = cat(1,x0,x);
y = cat(1,y0,y);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % lower circular segment at lower right
% cc = [initTarget.xmax initTarget.ymin] ;
% rr = tMax*velocity;
% th = 3*pi/2:ds:2*pi;
% xy = [ones(size(th'))*cc(1) ones(size(th'))*cc(2)] + rr * [cos(th') sin(th')];
% x0 = xy(:,1);
% y0 = xy(:,2);
% 
% x = cat(1,x0,x);
% y = cat(1,y0,y);

% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Top edge
% x0 = initTarget.xmin:ds:initTarget.xmax;
% y0 = (initTarget.ymax + vtarget*tMax)*ones(size(x0));
% x = cat(1,x0',x);
% y = cat(1,y0',y);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Left top slanted part
ptcy = initTarget.ymax + vtarget*tMax;
ptcx = initTarget.xmin;
m = sqrt(vtarget^2-velocity^2) / velocity;
xa = ptcx-vtarget*tMax/(m^-1+m); % coordinates of beginning of arc
ya = ptcy-vtarget*tMax/(m^-2+1);
x0 = xa:0.001:ptcx;
y0 = ya + m*(x0-xa);

% Refine resolution
plength = sum(sqrt(diff(x0).^2 + diff(y0).^2));   % length of segment
% ds = plength/n;                              % want 1000 points per 0.1 segment length
x0 = linspace(xa,ptcx,plength*n);
y0 = ya + m*(x0-xa);

x = cat(1,x0',x);
y = cat(1,y0',y);

thetaM = pi;
thetam = atan2(ya - initTarget.ymax, xa-ptcx);
theta = thetam:0.001:thetaM;

x0 = ptcx+ velocity*tMax*cos(theta);
y0 = initTarget.ymax + velocity*tMax*sin(theta);

% Refine resolution
plength = sum(sqrt(diff(x0).^2 + diff(y0).^2));   % length of segment
% ds = plength/n;                              % want 1000 points per 0.1 segment length

theta = linspace(thetam,thetaM,round(plength*n));

x0 = ptcx+ velocity*tMax*cos(theta);
y0 = initTarget.ymax + velocity*tMax*sin(theta);

x = cat(1,x0',x);
y = cat(1,y0',y);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Right top slanted part
% ptcx = initTarget.xmax;
% m = -sqrt(vtarget^2-velocity^2) / velocity;
% xa = ptcx-vtarget*tMax/(m^-1+m); % coordinates of beginning of arc
% ya = ptcy-vtarget*tMax/(m^-2+1);
% x0 = ptcx:ds:xa;
% y0 = ya + m*(x0-xa);
% x = cat(1,x0',x);
% y = cat(1,y0',y);
% 
% thetam = 0;
% thetaM = atan2(ya - initTarget.ymax, xa-ptcx);
% theta = thetam:ds:thetaM;
% 
% x0 = ptcx+ velocity*tMax*cos(theta);
% y0 = initTarget.ymax + velocity*tMax*sin(theta);
% x = cat(1,x0',x);
% y = cat(1,y0',y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obstacle shadow
ptcy = initObs.ymin + vobs*tMax;
ptcx = initObs.xmin;
m = -sqrt(vobs^2-velocity^2) / velocity;
x0 = ptcx:0.001:0;
y0 = ptcy + m*(x0-ptcx);

% Refine resolution
plength = sum(sqrt(diff(x0).^2 + diff(y0).^2));   % length of segment
% ds = plength/n;                              % want 1000 points per 0.1 segment length

x0 = linspace(ptcx,0,round(plength*n));
y0 = ptcy + m*(x0-ptcx);

x = cat(1,x0',x);
y = cat(1,y0',y);

% RIGHT SIDE
% ptcx = initObs.xmax;
% m = sqrt(vobs^2-velocity^2) / velocity;
% x0 = 0:ds:ptcx;
% y0 = ptcy + m*(x0-ptcx);
% x = cat(1,x0',x);
% y = cat(1,y0',y);

plot(x,y,'r.','markersize',1)
title(['t=' num2str(-tMax)])


e_avg = zeros(size(Nx));
e_max = zeros(size(Nx));
dx = zeros(size(Nx));

for i = 1:length(Nx)
    e = eval_u(g{i},sd_reach{i},[x y]);
    e_avg(i) = mean(e);
    e_max(i) = max(abs(e));
    dx(i) = mean(g{i}.dx);
end

figure; 
% plot(Nx,e_avg,'bo-'); hold on
% plot(Nx,e_max,'ro-'); hold on
% plot(Nx,dx,'k.-')
% semilogy(2./Nx,e_avg,'bo-'); hold on
% semilogy(2./Nx,e_max,'ro-'); hold on
% semilogy(2./Nx,dx,'k.-')
loglog(dx,e_avg,'bo-'); hold on
loglog(dx,e_max,'ro-'); hold on
loglog(dx,dx,'k--')
legend('e_{avg}','e_{max}','dx')
ylabel('Error')
xlabel('dx')
xlim([10^-2.2 10^-1])