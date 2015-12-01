% Simulates pursuit-evasion game

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
 x = [0.7; -0.8;2.5];

t_capture = inf;
freeze_at_capture = 0;
skipPlots = 0;

capture2D = sqrt((g2D.xs{1} - xEvader(x(3,end))).^2 + (g2D.xs{2} - yEvader(x(3,end))).^2) - captureRadius;

if ~skipPlots
    f1 = figure;
    f2 = figure;

    % Initial time

    i = 1;
    figure(f1),subplot(2,2,1) 
    plot (xEvader(g.vs{3}),yEvader(g.vs{3}),'--','Color',[.5,.5,.5]); hold on
    %contour(g2D.xs{1},g2D.xs{2},target2D(:,:,1),[0 0], 'r'); hold on
    contour(g2D.xs{1},g2D.xs{2},capture2D(:,:,1),[0 0], 'LineColor',[0,.75,0]);
    [~, value2D] = proj2D(g, [0 0 1], N2D, value(:,:,:,i), x(3,i));
    contour(g2D.xs{1},g2D.xs{2},value2D,[0 0], 'b');
    contour(g2D.xs{1},g2D.xs{2},obstacle2D(:,:,1),[0 0], 'k');
    plot(x(1,1),x(2,1),'^','MarkerEdgeColor','k','MarkerFaceColor','r');
    plot(xEvader(x(3,1)), yEvader(x(3,1)), 'o','MarkerEdgeColor','k','MarkerFaceColor','b')

    axis square; axis(g.axis(1:4))

    xvalue = eval_u(g,value(:,:,:,1),x);
    title(['v(x, ' num2str(tau(1)) ') = ' num2str(xvalue)])

    % Time progression

    figure(f1),subplot(2,2,2)
    plot (xEvader(g.vs{3}),yEvader(g.vs{3}),'--','Color',[.5,.5,.5]); hold on
    hp_now   = plot(x(1,end),x(2,end),'^','MarkerEdgeColor','k','MarkerFaceColor','r');
    hp_trace = plot(x(1,:),x(2,:),'r-');
    he_now   = plot(xEvader(x(3,end)), yEvader(x(3,end)), 'o','MarkerEdgeColor','k','MarkerFaceColor','b');
    he_trace = plot(xEvader(x(3,:)), yEvader(x(3,:)), 'b-');
    axis(g.axis(1:4))
    axis square

    obs_value = eval_u(g2D,obstacle2D(:,:,1),x(1:2,1));
    %tar_value = eval_u(g2D,target2D(:,:,1),x(:,:,1));
    tar_value = eval_u(g2D,capture2D(:,:,1),x(:,:,1));

    %[~, ht] = contour(g2D.xs{1},g2D.xs{2},target2D(:,:,1),[0 0], 'r');
    [~, ht] = contour(g2D.xs{1},g2D.xs{2},capture2D(:,:,1),[0 0], 'LineColor',[0,.75,0]);

    [~, ho] = contour(g2D.xs{1},g2D.xs{2},obstacle2D(:,:,1),[0 0], 'k');
    [~, value2D_i] = proj2D(g, [0 0 1], N2D, value(:,:,:,i), x(3,i));
    [~, hr] = contour(g2D.xs{1},g2D.xs{2},value2D(:,:,i),[0 0], 'b');

    subplot(2,2,3) % Plot target values
    hl = plot(tau(1:i), tar_value); hold on
    plot(tau, ones(size(tau))*xvalue,'k:')
    title('l(x(t), t)')
    xlim([tau(1) tau(end)])

    subplot(2,2,4) % Plot obstacle values
    hg = plot(tau(1:i),-obs_value); hold on
    plot(tau, zeros(size(tau)), 'k:')
    xlim([tau(1) tau(end)])
    title('g(x(t),t)')

    % Save snapshots as we go: code courtesy of Mo, Qie and Casey -------------
    tplot = [0, 0.3 0.55263, 1];
    numPlots = length(tplot);
    spC = ceil(sqrt(numPlots));
    spR = ceil(numPlots/spC);
    plotnum = 0;
    
% Begin Make Video --------------------------------------------------------
    %     % Prepare the new file.
    vidObj = VideoWriter('1stGameVideo','MPEG-4');
    vidObj.FrameRate = 5;
    vidObj.Quality = 100; % maximum quality
    open(vidObj);
    fv = figure;
    av = axes;
    axis(g.axis(1:4))
    axis square
    hp_now_v   = copyobj(hp_now,av);
    hp_now_v.MarkerSize = 16;
    hp_trace_v = copyobj(hp_trace,av);
    he_now_v   = copyobj(he_now,av);
    he_now_v.MarkerSize = 16;
    he_trace_v = copyobj(he_trace,av);
    ht_v = copyobj(ht,av);
    ho_v = copyobj(ho,av);
    hr_v = copyobj(hr,av);
    
 
%     Create an animation.
%     Z = peaks; surf(Z);
%     axis tight
%     set(gca,'nextplot','replacechildren');
%  
%     for k = 1:20
%        surf(sin(2*pi*k/20)*Z,Z)
%  
%        % Write each frame to the file.
%        currFrame = getframe;
%        writeVideo(vidObj,currFrame);
%     end
%   
%     % Close the file.
%     close(vidObj);
    
% End Make Video ----------------------------------------------------------  
    
end
%--------------------------------------------------------------------------
startTime = cputime;
for i = 2:length(tau)
    % Update state and save to trajectory
    P = extractCostates(g,value(:,:,:,i-1));
    p = calculateCostate(g,P,x(:,i-1),periodicDim);
    dir1 = -p(1:2)/norm(p(1:2));
    dir2 = p(3)/norm(p(3));
    if ~freeze_at_capture || t_capture==inf
        xdot = [(1-0.5*mean([tau(i-1),tau(i)]))*velocity1*dir1(:);velocity2*dir2(:)];
    else
        xdot = [0;0;0];
    end
    x = [x, x(:,i-1) + xdot*(tau(i)-tau(i-1))]; % Forward-Euler Integration
    
    if skipPlots,continue,end
    
    % Update obstacle and target values
    capture2D = cat(3,capture2D,...
        sqrt((g2D.xs{1} - xEvader(x(3,i))).^2 + (g2D.xs{2} - yEvader(x(3,i))).^2) - captureRadius);
    obs_value = [obs_value, eval_u(g2D,obstacle2D(:,:,i),x(:,i))];
    %tar_value = [tar_value, eval_u(g2D,target2D(:,:,i),x(:,i))];  
    tar_value = [tar_value, eval_u(g2D,capture2D(:,:,i),x(:,i))]; 
    
    if tar_value(i)<=0 && t_capture==inf
        t_capture = tau(i);
        disp(['Capture at t = ' num2str(tau(i))])
    end
    
    %%% Update Plot
    figure(f1),subplot(2,2,2)
    % Players
    hp_now.XData   = x(1,i);                hp_now.YData   = x(2,i);
    hp_trace.XData = x(1,:);                hp_trace.YData = x(2,:);
    he_now.XData   = xEvader(x(3,i));       he_now.YData   = yEvader(x(3,i));
    he_trace.XData = xEvader(x(3,:));       he_trace.YData = yEvader(x(3,:));
    
    % Level Sets
    [~, value2D_i] = proj2D(g, [0 0 1], N2D, value(:,:,:,i), x(3,i));
    value2D = cat(3,value2D,value2D_i);
    %ht.ZData = target2D(:,:,i);
    ht.ZData = capture2D(:,:,i);
    ho.ZData = obstacle2D(:,:,i);
    hr.ZData = value2D_i;
    
    % Title
    title(['l(x(' num2str(tau(i)) '), ' num2str(tau(i)) ') = ' num2str(tar_value(end))])
    
    subplot(2,2,3) % Plot target values
    hl.XData = tau(1:i);        hl.YData = tar_value;
    
    subplot(2,2,4) % Plot obstacle values
    hg.XData = tau(1:i);        hg.YData = -obs_value;
    
% Save snapshots as we go: code courtesy of Mo, Qie and Casey -------------
    if ~isempty(tplot) && tau(i) >= tplot(1)
        plotnum = plotnum+1;
        figure(f2)
        s = subplot(spR,spC,plotnum);
        copyobj(f1.Children(3).Children, s); % Child 3 is subplot (2,2,1)
        tplot(1) = [];
        
        title(sprintf('t=%0.2f', tau(i)));
        
        axis square; axis(g.axis(1:4))
        
    end
% Begin Make Video --------------------------------------------------------
    
    % Update video plot
    figure(fv)
    % Players
    hp_now_v.XData   = x(1,i);                hp_now_v.YData   = x(2,i);
    hp_trace_v.XData = x(1,:);                hp_trace_v.YData = x(2,:);
    he_now_v.XData   = xEvader(x(3,i));       he_now_v.YData   = yEvader(x(3,i));
    he_trace_v.XData = xEvader(x(3,:));       he_trace_v.YData = yEvader(x(3,:));
    
    % Level Sets
    ht_v.ZData = capture2D(:,:,i);
    ho_v.ZData = obstacle2D(:,:,i);
    hr_v.ZData = value2D_i;

    % Write each frame to the file.
    figure(fv)
    currFrame = getframe;
    writeVideo(vidObj,currFrame);
       
% End Make Video ----------------------------------------------------------

%--------------------------------------------------------------------------

end

% Begin Make Video --------------------------------------------------------
% Close the file.
     close(vidObj);
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
%% Create additional plots

% To be written tomorrow