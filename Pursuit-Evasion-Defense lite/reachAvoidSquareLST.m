function [ reach, g, target, tau, obstacle] = reachAvoidSquareLST(Nx, accuracy)
% air3D: demonstrate the 3D aircraft collision avoidance example
%
%   [ data, g, data0 ] = air3D(accuracy)
%
% In this example, the target set is a circle at the origin (cylinder in 3D)
% that represents a collision in relative coordinates between the evader
% (player a, fixed at the origin facing right) and the pursuer (player b).
%
% The relative coordinate dynamics are
%
%   \dot x    = -v_a + v_b \cos \psi + a y
%	  \dot y    = v_b \sin \psi - a x
%	  \dot \psi = b - a
%
% where v_a and v_b are constants, input a is trying to avoid the target
%	input b is trying to hit the target.
%
% For more details, see my PhD thesis, section 3.1.
%
% This function was originally designed as a script file, so most of the
% options can only be modified in the file.  For example, edit the file to
% change the grid dimension, boundary conditions, aircraft parameters, etc.
%
% To get exactly the result from the thesis choose:
%   targetRadius = 5, velocityA = velocityB = 5, inputA = inputB = +1.
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
tMax = 0.5;                  % End time.
plotSteps = 12;               % How many intermediate plots to produce?
t0 = 0;                      % Start time.
singleStep = 1;              % Plot at each timestep (overrides tPlot).

% Period at which intermediate plots should be produced.
tPlot = (tMax - t0) / (plotSteps - 1);

% How close (relative) do we need to get to tMax to be considered finished?
small = 100 * eps;

% What kind of dissipation?
dissType = 'global';

%---------------------------------------------------------------------------
% Problem Parameters.
velocity = 0.5;   % velocity of vehicle
restricted = 0; % If set to 1, vehicle can only choose directions between -pi/4 and pi/4
vtarget = 1.5;    % velocity of the target set (moving downwards in forward time)
vobs = 1;       % velocity of obstacle (moving left in forward time)


% Pause after each plot?
pauseAfterPlot = 0;

% Delete previous plot before showing next?
deleteLastPlot = 1;

% Plot in separate subplots (set deleteLastPlot = 0 in this case)?
useSubplots = 0;



%---------------------------------------------------------------------------
% Approximately how many grid cells?
%   (Slightly different grid cell counts will be chosen for each dimension.)
if (nargin < 1)
    Nx = 201;
end

if(nargin < 2)
    accuracy = 'veryHigh';
end



% Create the grid.
g.dim = 2;
g.min = [  -1; -1];
g.max = [ +1; +1 ];
g.bdry = { @addGhostExtrapolate; @addGhostExtrapolate};
% Roughly equal dx in x and y (so different N).
g.N = [ Nx; Nx ];
% Need to trim max bound in \psi (since the BC are periodic in this dimension).
g = processGrid(g);

%---------------------------------------------------------------------------
% Create initial conditions (cylinder centered on origin).
% data = shapeRectangleByCorners(g,[0 -0.25+tMax], [0.3 0.05+tMax]);
initTarget.xmin = -0.2;
initTarget.xmax = 0.2;
initTarget.ymin = -0.2;
initTarget.ymax = 0.2;

data = shapeRectangleByCorners(g,[initTarget.xmin initTarget.ymin], [initTarget.xmax initTarget.ymax]);
target = data;


initObs.xmin = -0.1;
initObs.xmax = 0.1;
initObs.ymin = -0.6;
initObs.ymax = -0.4;
obstacle = shapeRectangleByCorners(g,[initObs.xmin initObs.ymin], [initObs.xmax initObs.ymax]);
tau = t0;

reach = shapeDifference(target,obstacle);
data = shapeDifference(target,obstacle);
%---------------------------------------------------------------------------
% Set up spatial approximation scheme.
schemeFunc = @termLaxFriedrichs;
schemeData.hamFunc = @RASHamFunc;
schemeData.partialFunc = @RASPartialFunc;
schemeData.grid = g;

% The Hamiltonian and partial functions need problem parameters.
schemeData.velocity = velocity;
schemeData.restricted = restricted;
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
integratorOptions = odeCFLset('factorCFL', 0.25, 'stats', 'off');

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

h = cell(1);
[~, h{1}] = contour(g.xs{1},g.xs{2}, data, [0 0], 'b');

hold on;
axis(g.axis);
axis square

drawnow;

%---------------------------------------------------------------------------
% Loop until tMax (subject to a little roundoff).
tNow = t0;
startTime = cputime;
k = 0;
while(tMax - tNow > small * tMax)
    k = k+1;
    % Reshape data array into column vector for ode solver call.
    y0 = data(:);
    
    % How far to step?
    tNext = min(tMax, tNow + tPlot);
    tSpan = [ tNow, tNext ];

    %------------------------------------------------------------------------
    % Set up masking so that the reachable set does not propagate through the
    % avoid set
%     movingTarget = shapeRectangleByCorners(g,[0 -0.25+tMax-tNext], [0.3 0.05+tMax-tNext]);
%     schemeData.maskData = movingTargetV;
%     schemeData.maskFunc = @min;
    
    % Take a timestep.
    [ t y ] = feval(integratorFunc, schemeFunc, tSpan, y0,...
        integratorOptions, schemeData);
    tNow = t(end);

    % Moving obstacle
    
%     movingObsV = movingObs(:);
%     y = min(y,movingTargetV);
%     y(movingObsV <= 0) = 1e5;
    
    tau = cat(1,-tNow,tau);
    % Get back the correctly shaped data array
    data = reshape(y, g.shape);

    movingTarget = shapeRectangleByCorners(g,[initTarget.xmin initTarget.ymin+vtarget*tNow], [initTarget.xmax initTarget.ymax+vtarget*tNow]);
    target = cat(3,movingTarget,target);
    
    movingObs = shapeRectangleByCorners(g,[initObs.xmin initObs.ymin+vobs*tNow], [initObs.xmax initObs.ymax+vobs*tNow]);
    obstacle = cat(3,movingObs,obstacle);
%     movingObs = shapeRectangleByCorners(g,[-0.3+vobs*tNow 0], [-0.2+vobs*tNow 0.1]);
%     movingObs = ones(g.shape);
%     movingObs = shapeRectangleByCorners(g,[-0.3 -0.2], [-0.15 0.1]);
    
    data = min(data, movingTarget);
    data = max(data, -movingObs);
    
    
    reach = cat(3,data,reach);
    if(pauseAfterPlot)
        % Wait for last plot to be digested.
        pause;
    end
    
    % Get correct figure, and remember its current view.
    figure(f);
    
    % Delete last visualization if necessary.
    if(deleteLastPlot)
        for i = 1:length(h)
            delete(h{i});
        end
    end
    
    % Move to next subplot if necessary.
    if(useSubplots)
        plotNum = plotNum + 1;
        subplot(rows, cols, plotNum);
    end
    
    % Create new visualization.
%     h = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(tNow) ]);
    [~, h{1}] = contour(g.xs{1},g.xs{2}, data, [0 0], 'b');
%     [~, h{1}] = contour(g.xs{1},g.xs{2}, data, 0:0.1:1);
    hold on
    [~, h{2}] = contour(g.xs{1},g.xs{2},movingTarget,[0 0], 'r');
    [~, h{3}] = contour(g.xs{1},g.xs{2},movingObs,[0 0], 'k');
    title(['t=' num2str(tNow, '%4.3f')])
    axis square
    
    drawnow;
    export_fig(['2Dsim/' num2str(k)], '-png')
end

endTime = cputime;
fprintf('Total execution time %g seconds\n', endTime - startTime);


%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function hamValue = RASHamFunc(t, data, deriv, schemeData)

checkStructureFields(schemeData, 'grid', 'velocity', 'restricted');

g = schemeData.grid;
v = schemeData.velocity;


abs_Du = sqrt(deriv{1}.^2 + deriv{2}.^2);

if schemeData.restricted    % Only allowed to move in -pi/4 to pi/4 directions
    ind1 = logical( deriv{1}<=0 & abs(deriv{2})<=abs(deriv{1}) );   
    ind2 = ~ind1;
    hamValue = zeros(size(data));
    hamValue(ind1) = abs_Du(ind1);
    hamValue(ind2) = -1/sqrt(2) * deriv{1}(ind2) + 1/sqrt(2) * abs(deriv{2}(ind2));
    hamValue = v*hamValue;
    hamValue = max(0, hamValue);
else                        % allowed to move in any direction
    hamValue = v*sqrt(deriv{1}.^2 + deriv{2}.^2);
end

% hamValue = -hamValue;





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

checkStructureFields(schemeData, 'grid', 'velocity');

v = schemeData.velocity;

norm_v = sqrt(derivMax{1}.^2 + derivMax{2}.^2);

switch dim
    case 1
        alpha = v * abs(derivMax{1} ./ norm_v);
        
    case 2
        alpha = v * abs(derivMax{2} ./ norm_v);

end
