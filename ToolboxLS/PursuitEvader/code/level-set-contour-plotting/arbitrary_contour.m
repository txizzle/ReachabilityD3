% Plot contours of some arbitrary function z(x,y)
% Feel free to play around and modify z to get some intuition

% Define "grid"
xx=0:0.1:5;
yy=0:0.1:3;

% Define function
zz = zeros(length(xx),length(yy));
for i=1:length(xx)
    for j=1:length(yy)
        zz(i,j) = yy(j)^3 - xx(i)^2; % Define z(x,y) to be whatever function
    end
end

% Plot contours
figure;
contour_values = [0,1,2,3]; % If you only want one value, include it twice, e.g. [0,0]
[C,~] = contour(xx,yy,zz',contour_values);
xlabel('x','FontSize',14),ylabel('y','FontSize',14)
title('Contour plot for z(x,y)','FontSize',12,'FontWeight','normal');
axis equal


%% Notes

% The convention MATLAB uses for this function ("help contour") is:
% - the xx vector contains the horizontal coordinates of the grid
% - the yy vector contains the horizontal coordinates of the grid
% - the zz matrix has x as the column index and y as the row index, i.e.
%
%               x1  x2  x3  x4  x5
%           y1   *   *   *   *   *
%   zz1 =   y2   *   *   *   *   *
%           y3   *   *   *   *   *
%
% This is a slight problem for us because our standard convention is that
% the zz (hyper)matrix has first x, then y and so on for higher dimensions.
% That, is we represent our data as:
%
%               y1  y2  y3
%           x1   *   *   *
%           x2   *   *   *
%   zz2 =   x3   *   *   *
%           x4   *   *   *
%           x5   *   *   *
%
% Eventually, the picture needs to have the x increasing left to right and
% the y increasing bottom to top, that is:
%
%           y3   *   *   *   *   *
%   pic =   y2   *   *   *   *   *
%           y1   *   *   *   *   *
%               x1  x2  x3  x4  x5
%
% Since malab produces a contour plot in the form "pic" assuming "zz1",
% all we need to do is feed it our "zz2" transposed, as in the code above.


% Also, from "help contour":
%
% [C,H] = contour(...) returns contour matrix C and a handle, H, to
% a contour object.
%
% ...which means we can feed this (lighter) information to javascript.