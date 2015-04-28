function xs = xEvader(ss)
% This function computes the horizontal position of the evader from the
% curve parametrization s of a rectangle defined as:
% Bottom side:   0<s<1.2
% Right  side: 1.2<s<2.8
% Top    side: 2.8<s<4.0
% Left   side: 4.0<s<5.6(=0)

xs = arrayfun(@xScalarEvader,ss);

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