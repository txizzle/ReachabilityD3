function xy = xyEvader(xs)
% This function computes the position of the evader in (x,y) space from the
% curve parametrization s of a rectangle defined as:
% Bottom side:   0<s<1.2
% Right  side: 1.2<s<2.8
% Top    side: 2.8<s<4.0
% Left   side: 4.0<s<5.6(=0)

xy{1} = arrayfun(@xScalarEvader,xs);
xy{2} = arrayfun(@yScalarEvader,xs);

end


% Auxiliary functions:

function x = xScalarEvader(s)

% Function takes in each grid point

if s<1.2
    x = s-0.6;
elseif s<2.8
    x = 0.6;
elseif s<4.0
    x = 0.6 - (s-2.8);
else
    x = -0.6;
end

end


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