function ys = yEvader(ss)
% This function computes the vertical position of the evader from the
% curve parametrization s of a rectangle defined as:
% Bottom side:   0<s<1.2
% Right  side: 1.2<s<2.8
% Top    side: 2.8<s<4.0
% Left   side: 4.0<s<5.6(=0)

ys = arrayfun(@yScalarEvader,ss);

end

% Auxiliary functions:

function y = yScalarEvader(s)

% Function takes in each grid point

if s<1.2
    y = -0.8;
elseif s<2.8
    y = -0.8 + (s-1.2);
elseif s<4.0
    y = 0.8;
else
    y = 0.8 - (s-4.0);
end

end