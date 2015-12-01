clear all; close all;
mex movingTargetReachAvoid4D.cpp

%% GRID SETUP
Nx = 45;

% Create the computation grid.
g.dim = 4;
g.min = [  -1; -1; -1; -1];
g.max = [ +1; +1; +1; +1 ];
g.bdry = { @addGhostExtrapolate; @addGhostExtrapolate; @addGhostExtrapolate; @addGhostExtrapolate };
g.N = [ Nx; Nx; Nx; Nx ];
g = processGrid(g);

%% GAME SETUP
game = 'midTarget_LObs_vJ';
% game = 'OLGameModified';
run(game);
load([game '_4DHJI'])

% Initial conditions
u = shapeUnion(attackerWin,obs_d);
u0 = shapeUnion(attackerWin,obs_d);

% Avoid set
a = shapeUnion(defenderWin,obs_a);

h = g.dx(1);
k = 0.005;
tmax = 3;
t = k:k:tmax;
T = length(t);

tic; movingTargetReachAvoid4D(u, a, h, k, T, velocitya, velocityd);  toc;

% Fix attacker
figure
for i = 1:length(xa_init)
    subplot(2,2,i)
    hold off
    [~, hs] = visualizeGame(g2D, target2D, obs2D, xa_init{i}, xd_init, captureRadius, dom_map);
    
    [~, data2D{i}] = proj2D(g,[1 1 0 0], N2D, u, xa_init{i});
    [~, hua] = contour(g2D.xs{1}, g2D.xs{2},  data2D{i}, [0 0], 'r');
    
    [~, data2D{i}] = proj2D(g,[1 1 0 0], N2D, data, xa_init{i});
    [~, hdataa] = contour(g2D.xs{1}, g2D.xs{2},  data2D{i}, [0 0], 'b');
        
end
axis square
legend([hua hdataa], {'C++','MatLab'})

% Fix defender
figure
for i = 1:length(xd_init)
    subplot(2,2,i)
    hold off
    [~, hs] = visualizeGame(g2D, target2D, obs2D, xa_init, xd_init{i}, captureRadius, dom_map);
    
    [~, data2D{i}] = proj2D(g,[0 0 1 1], N2D, u, xd_init{i});
    [~, hud] = contour(g2D.xs{1}, g2D.xs{2},  data2D{i}, [0 0], 'r');
    
    [~, data2D{i}] = proj2D(g,[0 0 1 1], N2D, data, xd_init{i});
    [~, hdatad] = contour(g2D.xs{1}, g2D.xs{2},  data2D{i}, [0 0], 'b');    
end
axis square
legend([hud hdatad], {'C++','MatLab'})