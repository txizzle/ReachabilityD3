clear all; close all;
mex movingTargetReachAvoid.cpp

%% GRID SETUP
Nx = 201;

% Create the grid.
g.dim = 2;
g.min = [  -1; -1];
g.max = [ +1; +1];
g.bdry = { @addGhostExtrapolate; @addGhostExtrapolate};
% Roughly equal dx in x and y (so different N).
g.N = [ Nx; Nx ];
% Need to trim max bound in \psi (since the BC are periodic in this dimension).
g = processGrid(g);

%% INITIALIZATION
k = 0.005;
T = 0:k:0.5;
[~,~,t] = ndgrid(g.vs{1},g.vs{2},T);

% Initialize moving target and obstacle
u1 = zeros(size(t));
l = zeros(size(t));
% a = zeros(size(t));
a = ones(size(t));
for i = 1:length(T)
    l(:,:,i) = shapeRectangleByCorners(g,[0 -0.25+T(i)], [0.3 0.05+T(i)]);
    a(:,:,i) = shapeRectangleByCorners(g,[-0.3+1*T(i) -0.2], [-0.15+1*T(i) -0.0]);
end

% return
A = 0.5;

h = double(g.dx(1));

%% SOLVE
% Solve moving target and obstacle
u1(:,:,1) = l(:,:,1);
tic; movingTargetReachAvoid(u1, l, a, h, k, A); toc

%% PLOT
% Plot moving target reachable set (for all available time horizons)
numPlots = 12;
subC = ceil(sqrt(numPlots));
subR = ceil(numPlots/subC);
pInd = round(linspace(1, length(T), numPlots));

h1 = figure;
for i = 1:length(pInd)
    subplot(subR, subC, i)
    contour(g.xs{1},g.xs{2},u1(:,:,pInd(i)), [0 0], 'linecolor','r'); hold on
    contour(g.xs{1},g.xs{2},a(:,:,pInd(i)), [0 0], 'linecolor','k'); hold on
    contour(g.xs{1},g.xs{2},l(:,:,pInd(i)), [0 0], 'linecolor','b'); hold off
       
    title(['t=' num2str(T(pInd(i)))])
%     pause(0.1)
%     pause
    axis square
    drawnow;
    grid on
end