clear all

[ reach, g, target, tau, obstacle] = reachAvoidSquareLST;


velocity = 0.5;

xmin = -0.5;
xmax = 0.5;
ymin = -0.6;
ymax = 0.6;
N = 10;

% x = [-0.3 0];
% x = [xmin+(xmax-xmin)*rand(N,1) ymin+(ymax-ymin)*rand(N,1)];
% x = [-0.2 -0.4];
x = [-0.13 -0.43];

figure;
subplot(2,2,1) % Plot reachable set at final time
contour(g.xs{1},g.xs{2},target(:,:,1),[0 0], 'r'); hold on
contour(g.xs{1},g.xs{2},reach(:,:,1),[0 0], 'b');
contour(g.xs{1},g.xs{2},obstacle(:,:,1),[0 0], 'k');
plot(x(:,1),x(:,2),'b.');

axis square; axis(g.axis)

value = eval_u(g,reach(:,:,1),x);
if size(x,1) == 1
    title(['v(x, ' num2str(tau(1)) ') = ' num2str(value)])
end

subplot(2,2,2)
hv = plot(x(:,1),x(:,2),'b.'); hold on
ht = plot(x(:,1),x(:,2),'b-'); hold on

obs_value = eval_u(g,obstacle(:,:,1),x(:,:,1));
tar_value = eval_u(g,target(:,:,1),x(:,:,1));

[~, hmt] = contour(g.xs{1},g.xs{2},target(:,:,1),[0 0], 'r');
[~, hmo] = contour(g.xs{1},g.xs{2},obstacle(:,:,1),[0 0], 'k');
[~, hmr] = contour(g.xs{1},g.xs{2},reach(:,:,1),[0 0], 'b');
    
for i = 2:length(tau)
    % Update state and save to trajectory
    P = extractCostates(g,reach(:,:,i));
    p = calculateCostate(g,P,x(:,:,1));
    dir = -p/norm(p);
    x = cat(3, x(:,:,1) + velocity*dir*(tau(i)-tau(i-1)), x);
    
    % Update obstacle and target values
    obs_value = cat(1,obs_value,eval_u(g,obstacle(:,:,i),x(:,:,1)));
    tar_value = cat(1,tar_value,eval_u(g,target(:,:,i),x(:,:,1)));  
    
    % Plot
    subplot(2,2,2)
    
    delete(hv);     hv = plot(x(:,1,1),x(:,2,1),'b.');
    delete(hmt);    [~, hmt] = contour(g.xs{1},g.xs{2},target(:,:,i),[0 0], 'r');
    delete(hmo);    [~, hmo] = contour(g.xs{1},g.xs{2},obstacle(:,:,i),[0 0], 'k');
    delete(hmr);    [~, hmr] = contour(g.xs{1},g.xs{2},reach(:,:,i),[0 0], 'b');
    delete(ht);     ht = plot(squeeze(x(:,1,:)),squeeze(x(:,2,:)),'b:');
    
    if size(x,1)==1
        title(['l(x(' num2str(tau(i)) '), ' num2str(tau(i)) ') = ' num2str(tar_value(end))])
    end
    axis(g.axis)
    axis square
    drawnow;
    
    subplot(2,2,3) % Plot target values
    plot(tau(1:i), tar_value); hold on
    plot(tau, ones(size(tau))*value,'k:')
    title('l(x(t), t)')
    xlim([tau(1) tau(end)])
    
    subplot(2,2,4) % Plot obstacle values
    plot(tau(1:i),obs_value); hold on
    plot(tau, zeros(size(tau)), 'k:')
    xlim([tau(1) tau(end)])
    title('g(x(t),t)')
end

% Path length
dx = diff(x(1,1,:));
dy = diff(x(1,2,:));
path_length = sum(sqrt(dx.^2 + dy.^2));
disp(['Path length = ' num2str(path_length)])