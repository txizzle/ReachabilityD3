clear all; close all;
mex staticTargetReach.cpp
mex movingTargetReach.cpp

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

% Initialize moving obstacle
u1 = zeros(size(t));
l = zeros(size(t));
for i = 1:length(T)
    l(:,:,i) = shapeRectangleByCorners(g,[0 -0.25+T(i)], [0.3 0.05+T(i)]);
end

% Initialize stationary obstacles
u2 = cell(length(T),1);
for i = 1:length(T)
    u2{i} = zeros(size(t,1), size(t,2), size(t,3)-i+1);
    u2{i}(:,:,1) = l(:,:,i);
end

% return
A = 0.5;
epsilon = 0;

h = double(g.dx(1));
% return

%% SOLVE
% Solve moving target
u1(:,:,1) = l(:,:,1);
tic; movingTargetReach(u1, l, h, k, A, epsilon); toc

% % Solve static target using various initial conditions
% u3 = u2{end};
% for i = 1:length(T)-1
%     tic; staticTargetReach(u2{i}, h, k, A, epsilon); toc
%     u3 = shapeUnion(u3, u2{i}(:,:,end));
% end

%% PLOT
% Plot moving target reachable set (for all available time horizons)
h1 = figure;
for i = 1:length(T)
    figure(h1);
    contour(g.xs{1},g.xs{2},u1(:,:,i), [0 0], 'linecolor','r'); hold on
    contour(g.xs{1},g.xs{2},l(:,:,i), [0 0], 'linecolor','b'); hold off
       
    title(['t=' num2str(T(i))])
%     pause(0.1)
%     pause
    drawnow;
end

% % Plot union of static target reachable sets
% hold on;
% contour(g.xs{1}, g.xs{2}, u3, [0 0], 'linecolor','k')
% contour(x,y,u(:,:,end),[-0.099 -0.099])