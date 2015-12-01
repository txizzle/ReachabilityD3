function UTM_Safe_Flight()

% Program: Play Pursuit-Evasion
% Author:  Jaime F. Fisac
% Version: 1.0
% Date:    2015/07/26
%
% This code allows a user to control an aircraft in a reach-avoid scenario
% by entering the aircraft control input through the keyboard.
%
% If the option HJI_override is set to TRUE then the human input will
% be overridden by the optimal action to ensure mission success when
% necessary. If the option HJI_autonomous is set to TRUE then the pursuer
% always uses the policy from the HJI computation, ignoring user input.
%
% User input is processed through the figure object and therefore the
% corresponding information is stored and retreived using the figure
% itself, through guidata(). Any implemented code using this interface
% will need to keep this structure.
%
% USER CONTROLS: use either 'w','a','s','d' or the arrow keys
% while the figure is in focus to control the motion of the TRIANGLE on the
% screen. Use keys 'x' or 'esc' to close the window and end the simulation.


clear all
close all

solution   = load('1stGame_31');

tau        = solution.time_trace;
g          = solution.g;
g2D        = solution.g2D;
N2D        = size(g2D.xs{1})';
target2D   = solution.target_trace;
obstacle2D = solution.obstacle_trace;
value      = solution.value_trace;

periodicDim = [0 0 1]; % The evader moves along a closed path

captureRadius = 0.1;

%%% IMPORTANT: dynamics must be consistent with those used in pursuitEvasion3D
velocity1 = 3;
velocity2 = 2;

xmin = -0.5;
xmax = 0.5;
ymin = -0.6;
ymax = 0.6;
N    = 10;

% state is vertical, time trace grows horizontally
% x = [-0.3 0];
% x = [xmin+(xmax-xmin)*rand(N,1) ymin+(ymax-ymin)*rand(N,1)];
% x = [-0.2 -0.4];
% x = [-0.5; 0; 3];  % Easy capture for pursuer (interior of winning set)
% x = [-.8;-.25;3]; % Close escape (ish)
% x = [-.8;-.1;3];  % Capture at last instant
% x = [0.7; -0.8;2.5]; % Escape by fooling pursuer into obstacle
 x = [0.7; -0.4;2.5]; % Easy capture for pursuer from right side

t_capture = inf;
freeze_at_capture = 0;
HJI_override = 1;
HJI_autonomous = 0;

capture2D = sqrt((g2D.xs{1} - xEvader(x(3,end))).^2 + (g2D.xs{2} - yEvader(x(3,end))).^2) - captureRadius;



% Initial time

i = 1;
[~, value2D] = proj2D(g, [0 0 1], N2D, value(:,:,:,i), x(3,i));

%% Graphical setup
vidObj = VideoWriter('1stGameVideo','MPEG-4');
vidObj.FrameRate = 5;
vidObj.Quality = 100; % maximum quality
open(vidObj);
S.fh = figure('units','pixels',...
          ...'position',[500 500 300 300],...
          'menubar','none',...
          'name','pursuit-evasion',...
          'numbertitle','off',...
          'keypressfcn',@cmd_press,...   % callback for key press
          'keyreleasefcn',@cmd_release); % callback for key release    
av = axes;
axis(g.axis(1:4))
axis square
hold on
hp_now   = plot(x(1,end),x(2,end),'^','MarkerEdgeColor','k','MarkerFaceColor','r');
hp_trace = plot(x(1,:),x(2,:),'r-');
he_now   = plot(xEvader(x(3,end)), yEvader(x(3,end)), 'o','MarkerEdgeColor','k','MarkerFaceColor','b');
he_trace = plot(xEvader(x(3,:)), yEvader(x(3,:)), 'b-');
he_now.MarkerSize = 14;
hp_now.MarkerSize = 14;
[~, value2D] = proj2D(g, [0 0 1], N2D, value(:,:,:,i), x(3,i));
[~, ht] = contour(g2D.xs{1},g2D.xs{2},capture2D(:,:,1),[0 0], 'LineColor',[0,.75,0]);
[~, ho] = contour(g2D.xs{1},g2D.xs{2},obstacle2D(:,:,1),[0 0], 'k');
[~, hr] = contour(g2D.xs{1},g2D.xs{2},value2D(:,:,i),[0 0], 'b');
guidata(S.fh,S)
    
%% Define command keys
N_cmds = 5; % Here we are using 5 different controls (including 'quit')
S.cmdStatus = false(N_cmds,1);

% Allow alternative sets of control keys (for convenience)
S.cmdNames = {  'w', 'uparrow';
                's', 'downarrow';
                'a', 'leftarrow';
                'd', 'rightarrow';
                'x', 'escape'};
S.cmd.up    = 1;
S.cmd.down  = 2;
S.cmd.left  = 3;
S.cmd.right = 4;
S.cmd.quit  = 5;

guidata(S.fh,S) % store GUI information in figure object

%--------------------------------------------------------------------------
%% Main simulation loop
startTime = cputime;
for i = 2:length(tau)
    % Retrieve updated info from figure object
    S = guidata(S.fh);
    if S.cmdStatus(S.cmd.quit)
        fprintf('\nSimulation terminated by user.\n')
        break;
    end
    % Determine safety level
    RA_value = eval_u(g2D,value2D(:,:,i-1),x(:,i-1));
    % Update state and save to trajectory
    P = extractCostates(g,value(:,:,:,i-1));
    p = calculateCostate(g,P,x(:,i-1),periodicDim);
    if (RA_value < -0.01 || ~HJI_override) && ~HJI_autonomous
        if S.cmdStatus(S.cmd.up) && ~S.cmdStatus(S.cmd.down)
            dir1(2) = 0.1;
        elseif S.cmdStatus(S.cmd.down) && ~S.cmdStatus(S.cmd.up)
            dir1(2) = -0.1;
        else
            dir1(2) = 0;
        end
        if S.cmdStatus(S.cmd.left) && ~S.cmdStatus(S.cmd.right)
            dir1(1) = -0.1;
        elseif S.cmdStatus(S.cmd.right) && ~S.cmdStatus(S.cmd.left)
            dir1(1) = 0.1;
        else
            dir1(1) = 0;
        end
        if norm(dir1)
            dir1 = dir1/norm(dir1); % normalize if nonzero
        end
    else
        dir1 = -p(1:2)/norm(p(1:2)); % optimal (HJI) input
    end
    dir2 = p(3)/norm(p(3));
    if ~freeze_at_capture || t_capture==inf
        xdot = [(1-0.5*mean([tau(i-1),tau(i)]))*velocity1*dir1(:);velocity2*dir2(:)];
    else
        xdot = [0;0;0];
    end
    x = [x, x(:,i-1) + xdot*(tau(i)-tau(i-1))]; % Forward-Euler Integration
    
    
    %if skipPlots,continue,end
    
    % Update obstacle and target values
     capture2D = cat(3,capture2D,...
         sqrt((g2D.xs{1} - xEvader(x(3,i))).^2 + (g2D.xs{2} - yEvader(x(3,i))).^2) - captureRadius);

     [~, value2D_i] = proj2D(g, [0 0 1], N2D, value(:,:,:,i), x(3,i));
     value2D = cat(3,value2D,value2D_i);
     
    tar_value = eval_u(g2D,capture2D(:,:,i),x(:,i));
    if tar_value <=0 && abs(tau(i) - t_capture) > 0.1
        t_capture = tau(i);
        disp(['Capture at t = ' num2str(tau(i))])
    end

% Begin Make Video --------------------------------------------------------
    
    % Update video plot
    figure(S.fh)
    % Players
    hp_now.XData   = x(1,i);                hp_now.YData   = x(2,i);
    hp_trace.XData = x(1,:);                hp_trace.YData = x(2,:);
    he_now.XData   = xEvader(x(3,i));       he_now.YData   = yEvader(x(3,i));
    he_trace.XData = xEvader(x(3,:));       he_trace.YData = yEvader(x(3,:));
    
    % Level Sets
    ht.ZData = capture2D(:,:,i);
    ho.ZData = obstacle2D(:,:,i);
    hr.ZData = value2D_i;

    % Write each frame to the file.
    figure(S.fh)
    currFrame = getframe;
    writeVideo(vidObj,currFrame);
       
% End Make Video ----------------------------------------------------------

%--------------------------------------------------------------------------

end

% Begin Make Video --------------------------------------------------------
close(vidObj); % Close the file.
% End Make Video ----------------------------------------------------------


% Path lengths
dx_p = diff(x(1,:));
dy_p = diff(x(2,:));
path_length_p = sum(sqrt(dx_p.^2 + dy_p.^2));
ds_e = diff(x(3,:));
path_length_e = sum(abs(ds_e));
disp(['Pursuer path length = ' num2str(path_length_p)]);
disp(['Evader path length = ' num2str(path_length_e)]);

endTime = cputime;
compTime = endTime-startTime;
disp(['Computation time = ' num2str(compTime)]);


end

%% Key press and release callbacks
% These functions update the control input whenever keyboard status changes


function cmd_press(hObject,event)   
% This function updates the cmdStatus vector indicating what controls are
% being commanded by the user. It is automatically called when a key is
% pressed down on the keyboard.

S = guidata(hObject); % retrieve updated info from figure object
E = event;

S.cmdStatus = (any( strcmp(E.Key, S.cmdNames),2 ) | S.cmdStatus);
guidata(S.fh,S); % store update info in figure object

end


function cmd_release(hObject,event)
% This function updates the cmdStatus vector indicating what controls are
% being commanded by the user. It is automatically called when a key is
% released on the keyboard.

S = guidata(hObject); % retrieve updated info from figure object
E = event;

S.cmdStatus = (~any( strcmp(E.Key, S.cmdNames),2 ) & S.cmdStatus);
guidata(S.fh,S); % store updated info in figure object

end
