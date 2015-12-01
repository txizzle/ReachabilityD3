mex initFunction.cpp
x_min = 0;
y_min = -0.2;

x_max = 0.4;
y_max = 0.2;
ds = 0.01;

xtop = x_min:ds:x_max-ds;
ytop = y_max*ones(size(xtop));

yright = y_max:-ds:y_min+ds;
xright = x_max*ones(size(yright));

xbot = x_max:-ds:x_min+ds;
ybot = y_min*ones(size(xbot));

yleft = y_min:ds:y_max-ds;
xleft = x_min*ones(size(yleft));

xin = [xtop xright xbot xleft];
yin = [ytop yright ybot yleft];

figure;plot(xin,yin,'.')

h = 0.01;
X = -1:h:1;
Y = -1:h:1;

[x,y] = ndgrid(X,Y);
l = zeros(size(x));

initFunction(xin', yin', x, y, l)

hold on
contour(x,y,l,-0.1:0.025:1)