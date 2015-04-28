function b = obstacleMotion(t)
% Function for obstacle motion parametrized by single time-varying scalar.

if t<0.2
    b = -2*t;
elseif t<0.6
    b = -0.4;
elseif t<0.8
    b = -2.8 + 4*t;
else
    b = 0.4;
end