function [ g, g2D, time_trace, value_trace, target_trace, obstacle_trace,compTime] = pursuitEvasion3D(Nx, accuracy)
% air3D: demonstrate the 3D aircraft collision avoidance example
%
%   [ data, g, data0 ] = air3D(accuracy)
%
% System coordinates
%
% The state is given by [p_x, p_y, e_s]
%
% where e_s parametrizes the progress of the evader along a 1.2x1.6
% rectangle starting at the bottom left corner and moving counter-clockwise
% (the units of e_s are made to be equal to length units)
%
% Therefore we have:
% Bottom side:   0<s<1.2
% Right  side: 1.2<s<2.8
% Top    side: 2.8<s<4.0
% Left   side: 4.0<s<5.6(=0)
%
% The function xyEvader.m computes the true position (e_x, e_y) on a grid.
% Functions xEvader.m and yEvader.m do the same for each coordinate.
%
% Input Parameters:
%
%   accuracy: Controls the order of approximations.
%     'low': Use odeCFL1 and upwindFirstFirst.
%     'medium': Use odeCFL2 and upwindFirstENO2 (default).
%     'high': Use odeCFL3 and upwindFirstENO3.
%     'veryHigh': Use odeCFL3 and upwindFirstWENO5.
%
% Output Parameters:
%
%   data: Implicit surface function at t_max.
%
%   g: Grid structure on which data was computed.
%
%   data0: Implicit surface function at t_0.

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing
%   agreement contained in the file LICENSE in the top directory of
%   the distribution.
%
% Ian Mitchell, 3/26/04
% Subversion tags for version control purposes.
% $Date: 2012-07-04 14:27:00 -0700 (Wed, 04 Jul 2012) $
% $Id: air3D.m 74 2012-07-04 21:27:00Z mitchell $

%---------------------------------------------------------------------------
% You will see many executable lines that are commented out.
%   These are included to show some of the options available; modify
%   the commenting to modify the behavior.

%---------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
% run('addPathToKernel');

%---------------------------------------------------------------------------
% Integration parameters.
tMax = 1;                    % End time.
plotSteps = 12;              % How many intermediate plots to produce?
t0 = 0;                      % Start time.
singleStep = 1;              % Plot at each timestep (overrides tPlot).

% Period at which intermediate plots should be produced.
tPlot = (tMax - t0) / (plotSteps - 1); %%% CAREFUL: right now target&constraint only get updated at new plots!
                                       %%% FOR NOW: solve this by leaving singleStep = 1
% How close (relative) do we need to get to tMax to be considered finished?
small = 100 * eps;

% What kind of dissipation?
dissType = 'global';

%---------------------------------------------------------------------------
% Problem Parameters.
velocity1 = 3;   % velocity of pursuer
velocity2 = 2;   % velocity of evader
restricted = 0;  % If set to 1, vehicle can only choose directions between -pi/4 and pi/4
%vtarget = 1.5;   % velocity of the target set (moving downwards in forward time)
%vobs = 0.5;      % velocity of obstacle (moving left in forward time)


% Pause after each plot?
pauseAfterPlot = 0;

% Delete previous plot before showing next?
%deleteLastPlot = 1; % Updating handles instead

% Plot in separate subplots (set deleteLastPlot = 0 in this case)?
useSubplots = 0;



%---------------------------------------------------------------------------
% Approximately how many grid cells?
%   (Slightly different grid cell counts will be chosen for each dimension.)
if (nargin < 1)
    Nx = 51;
end

if(nargin < 2)
    accuracy = 'veryHigh';
end



% Create the grid.
g.dim = 3;
g.min = [  -1; -1; 0];
g.max = [ +1; +1; 5.6];
g.bdry = { @addGhostExtrapolate; @addGhostExtrapolate; @addGhostPeriodic};
% Roughly equal dx in x and y (so different N).
g.N = [ Nx; Nx; floor(5.6/2*Nx) ]; % Compensate for length of evader path vs side of square to keep density
% Need to trim max bound in \psi (since the BC are periodic in this dimension).
g = processGrid(g);

% Points where we will want to plot the evader for the 2D slices 
se = [0 1 2 3];

N2D = 200; % Number of grid points for the 2D grid when plotting (uses interpolation of value function)
g2D = proj2D(g, [0 0 1], N2D);

%---------------------------------------------------------------------------

% Initialize time signature
time_trace = tMax-t0;  % we make tau be forward time as opposed to backward

% data = shapeRectangleByCorners(g,[0 -0.25+tMax], [0.3 0.05+tMax]);
% initTarget.xmin = -0.2;
% initTarget.xmax =  0.2;
% initTarget.ymin = -0.8;
% initTarget.ymax = -0.6;

% data = shapeRectangleByCorners(g,[initTarget.xmin initTarget.ymin+vtarget*tMax -inf], [initTarget.xmax initTarget.ymax+vtarget*tMax inf]);
% target = shapeRectangleByCorners(g2D,[initTarget.xmin initTarget.ymin+vtarget*tMax], [initTarget.xmax initTarget.ymax+vtarget*tMax]);
% target3D = data;

%%% Obstacle
b = obstacleMotion(tMax); % terminal obstacle position
Obs.xmin = -0.2;
Obs.xmax =  0.2;
Obs.ymin = -0.6 + b;
Obs.ymax =  0.6 + b;
obstacleForPursuer = shapeRectangleByCorners(g,[Obs.xmin Obs.ymin -inf], [Obs.xmax Obs.ymax inf]);
obstacleForEvader  = interpn(g.xs{1}(:,1,1),g.xs{2}(1,:,1),obstacleForPursuer(:,:,1),xEvader(g.xs{3}),yEvader(g.xs{3}));

%%% Capture conditions
captureRadius = 0.1;
capture = sqrt( (g.xs{1} - xEvader(g.xs{3})).^2 + (g.xs{2} - yEvader(g.xs{3})).^2 ) - captureRadius;

%%% Implicit functions: target, constraint and value (old notation: data = value, reach = value_trace)
target = min(capture,obstacleForEvader);
constraint = max(-3,min(-4*obstacleForPursuer,3)); % "tighten" g(x,t) so that interior trajectories are close to optimal
value  = max(target, constraint); % formerly data
value_trace = value;              % formerly reach

%-----------------------
% Plotty stuff for later
capture2D = cell(4,1);
for i = 1:length(se)
    capture2D{i} = sqrt((g2D.xs{1} - xEvader(se(i))).^2 + (g2D.xs{2} - yEvader(se(i))).^2) - captureRadius;
end
obstacle_trace = shapeRectangleByCorners(g2D,[Obs.xmin Obs.ymin], [Obs.xmax Obs.ymax]);
target_trace = [];

%---------------------------------------------------------------------------
% Set up spatial approximation scheme.
schemeFunc = @termLaxFriedrichs;
schemeData.hamFunc = @RASHamFunc;
schemeData.partialFunc = @RASPartialFunc;
schemeData.grid = g;

% The Hamiltonian and partial functions need problem parameters.
schemeData.velocity1 = velocity1;
schemeData.velocity2 = velocity2;
schemeData.restricted = restricted;
schemeData.tMax = tMax;
%---------------------------------------------------------------------------
% Choose degree of dissipation.

switch(dissType)
    case 'global'
        schemeData.dissFunc = @artificialDissipationGLF;
    case 'local'
        schemeData.dissFunc = @artificialDissipationLLF;
    case 'locallocal'
        schemeData.dissFunc = @artificialDissipationLLLF;
    otherwise
        error('Unknown dissipation function %s', dissFunc);
end

%---------------------------------------------------------------------------


% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', 0.75, 'stats', 'on');

% Choose approximations at appropriate level of accuracy.
switch(accuracy)
    case 'low'
        schemeData.derivFunc = @upwindFirstFirst;
        integratorFunc = @odeCFL1;
    case 'medium'
        schemeData.derivFunc = @upwindFirstENO2;
        integratorFunc = @odeCFL2;
    case 'high'
        schemeData.derivFunc = @upwindFirstENO3;
        integratorFunc = @odeCFL3;
    case 'veryHigh'
        schemeData.derivFunc = @upwindFirstWENO5;
        integratorFunc = @odeCFL3;
    otherwise
        error('Unknown accuracy level %s', accuracy);
end

if(singleStep)
    integratorOptions = odeCFLset(integratorOptions, 'singleStep', 'on');
end

%---------------------------------------------------------------------------
% Initialize Display
f = figure;

% Set up subplot parameters if necessary.
if(useSubplots)
    cols = ceil(sqrt(plotSteps));
    rows = ceil(plotSteps / cols);
    plotNum = 1;
    subplot(rows, cols, plotNum);
end

h = cell(4,1);
ho = cell(4,1);
ht = cell(4,1);
for i = 1:length(se)
    subplot(2,2,i)
    [~, value2D] = proj2D(g, [0 0 1], N2D, value, se(i)); hold on
    [~, h{i}] = contour(g2D.xs{1},g2D.xs{2}, value2D, [0 0], 'b');

    plot(xEvader(se(i)), yEvader(se(i)), 'o','MarkerEdgeColor','k','MarkerFaceColor','b')
    contour(g2D.xs{1},g2D.xs{2}, capture2D{i}, [0 0],'LineColor',[0,.75,0]);

    %[~, ht{i}] = contour(g2D.xs{1},g2D.xs{2},movingTarget2D_trace, [0 0], 'r');
    [~, ho{i}] = contour(g2D.xs{1},g2D.xs{2},obstacle_trace,[0 0], 'k');
        
    axis(g.axis);
    axis square    
end
drawnow;






%---------------------------------------------------------------------------
% Loop until tMax (subject to a little roundoff).
tNow = t0;
startTime = cputime;
while(tMax - tNow > small * tMax)
    % Reshape data array into column vector for ode solver call.
    y0 = value(:);
    
    % How far to step?
    tNext = min(tMax, tNow + tPlot);    %%% CAREFUL: right now target&constraint only get updated at new plots!
    tSpan = [ tNow, tNext ];            %%% FOR NOW: solve this by leaving singleStep = 1

    %------------------------------------------------------------------------
    % Set up masking so that the reachable set does not propagate through the
    % avoid set
%     movingTarget = shapeRectangleByCorners(g,[0 -0.25+tMax-tNext], [0.3 0.05+tMax-tNext]);
%     schemeData.maskData = movingTargetV;
%     schemeData.maskFunc = @min;
    
    % Take a timestep.
    [ t, y ] = feval(integratorFunc, schemeFunc, tSpan, y0,...
        integratorOptions, schemeData);
    tNow = t(end);

    % Moving obstacle
    
%     movingObsV = movingObs(:);
%     y = min(y,movingTargetV);
%     y(movingObsV <= 0) = 1e5;
    
    time_trace = cat(1,tMax-tNow,time_trace);
    % Get back the correctly shaped data array
    value = reshape(y, g.shape);

    
    
%     movingTarget2D = shapeRectangleByCorners(g2D,[initTarget.xmin initTarget.ymin+vtarget*(tMax-tNow)], [initTarget.xmax initTarget.ymax+vtarget*(tMax-tNow)]);
%     movingTarget3D = shapeRectangleByCorners(g,[initTarget.xmin initTarget.ymin+vtarget*(tMax-tNow) -inf], [initTarget.xmax initTarget.ymax+vtarget*(tMax-tNow) inf]);
%     movingTarget2D_trace = cat(3,movingTarget2D,movingTarget2D_trace);
    
    %%% Obstacle
    b = obstacleMotion(time_trace(1));
    Obs.xmin = -0.2;
    Obs.xmax =  0.2;
    Obs.ymin = -0.6 + b;
    Obs.ymax =  0.6 + b;
    obstacleForPursuer = shapeRectangleByCorners(g,[Obs.xmin Obs.ymin -inf], [Obs.xmax Obs.ymax inf]);
    obstacleForEvader  = interpn(g.xs{1}(:,:,1),g.xs{2}(:,:,1),obstacleForPursuer(:,:,1),xEvader(g.xs{3}),yEvader(g.xs{3}));

    %%% Capture conditions
    % Unchanged (time-invariant)
    
    %%% Implicit functions: moving obstacle causes changes in target and constraint
    target = min(capture,obstacleForEvader);
    constraint = -max(-3,min(4*obstacleForPursuer,3));% "tighten" g(x,t) so that interior trajectories are close to optimal

    %%% Variational Inequality: value saturation step
    value = min(value, target);
    value = max(value, constraint);
    value_trace = cat(4,value,value_trace);
    
    
    %%% Begin plotty stuff
    movingObs2D = shapeRectangleByCorners(g2D,[Obs.xmin Obs.ymin], [Obs.xmax Obs.ymax]);
    obstacle_trace = cat(3, movingObs2D,obstacle_trace);
    if(pauseAfterPlot)
        % Wait for last plot to be digested.
        pause;
    end
    
    % Get correct figure, and remember its current view.
    figure(f);
    
    % Delete last visualization if necessary. %%% Updating handles instead
%     if(deleteLastPlot)
%         for i = 1:length(h)
%             delete(h{i});
%             delete(ht{i});
%             delete(ho{i});
%         end
%     end
    
    % Move to next subplot if necessary.
    if(useSubplots)
        plotNum = plotNum + 1;
        subplot(rows, cols, plotNum);
    end
    
    % Create new visualization.
%     h = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(tNow) ]);
    for i = 1:length(se)
        subplot(2,2,i)
        [~, value2D] = proj2D(g, [0 0 1], N2D, value, se(i));
        h{i}.ZData = value2D;

        %ht{i}.ZData = movingTarget2D;
        ho{i}.ZData = movingObs2D;
        title(['t=' num2str(time_trace(1))])
        axis square
    
        drawnow;        
    end
%     [~, h{1}] = contour(g.xs{1},g.xs{2}, data, 0:0.1:1);
%     hold on
    

    

end % end time loop

endTime = cputime;
compTime = endTime-startTime;
fprintf('Total execution time %g seconds\n', endTime - startTime);

figure
visualizeLevelSet(g, value, 'surface', 0); camlight

% Save data in .mat file
save(['1stGame_' num2str(Nx)], 'g', 'g2D', 'time_trace', 'value_trace', 'target_trace', 'obstacle_trace', 'compTime','-v7.3');


%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function hamValue = RASHamFunc(t_backward, data, deriv, schemeData)

checkStructureFields(schemeData, 'grid', 'velocity1', 'velocity2', 'restricted');

t = schemeData.tMax - t_backward;
g = schemeData.grid;
v1 = schemeData.velocity1;
v2 = schemeData.velocity2;

%hamValue = -1.5*t_backward*v1*sqrt(deriv{1}.^2 + deriv{2}.^2) + v2*abs(deriv{3});

% Time-invariant dynamics
%hamValue = -v1*sqrt(deriv{1}.^2 + deriv{2}.^2) + v2*abs(deriv{3});

% Pursuer slows down to half of its original speed
hamValue = -(1-0.5*t)*v1*sqrt(deriv{1}.^2 + deriv{2}.^2) + v2*abs(deriv{3});
hamValue = -hamValue;



%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function alpha = RASPartialFunc(t, data, derivMin, derivMax, schemeData, dim)
% air3DPartialFunc: Hamiltonian partial fcn for 3D collision avoidance example.
%
% alpha = air3DPartialFunc(t, data, derivMin, derivMax, schemeData, dim)
%
% This function implements the partialFunc prototype for the three dimensional
%   aircraft collision avoidance example (also called the game of
%   two identical vehicles).
%
% It calculates the extrema of the absolute value of the partials of the
%   analytic Hamiltonian with respect to the costate (gradient).
%
% Parameters:
%   t            Time at beginning of timestep (ignored).
%   data         Data array.
%   derivMin	 Cell vector of minimum values of the costate (\grad \phi).
%   derivMax	 Cell vector of maximum values of the costate (\grad \phi).
%   schemeData	 A structure (see below).
%   dim          Dimension in which the partial derivatives is taken.
%
%   alpha	 Maximum absolute value of the partial of the Hamiltonian
%		   with respect to the costate in dimension dim for the
%                  specified range of costate values (O&F equation 5.12).
%		   Note that alpha can (and should) be evaluated separately
%		   at each node of the grid.
%
%
% Ian Mitchell 3/26/04

checkStructureFields(schemeData, 'grid', 'velocity1', 'velocity2');

v1 = schemeData.velocity1;
v2 = schemeData.velocity2;

norm_v1 = sqrt(derivMax{1}.^2 + derivMax{2}.^2);

switch dim
    case 1
        alpha = v1 * abs(derivMax{1} ./ norm_v1);
        
    case 2
        alpha = v1 * abs(derivMax{2} ./ norm_v1);
        
    case 3
        alpha = v2;
end
