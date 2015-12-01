function v = eval_u(g, data, x)
% v = eval_u(g,data,x)
% Checks weather the point x is inside the set described by data
% 
% Inputs:
%   g       - grid
%   data    - implicit function describing the set
%   x       - points to check; each row is a point
%
% OUTPUT
%   v:  value at points x
%
% Mo Chen, 2014-09-16
% 
%%%----- modified by Jaime F. Fisac, 2015-03-23
%
% -- The switch is now over ndims(data) rather than g.dim
%
% -- If only one sample point is used, it can be in the form
% [x1; x2; ... ; xn] without a first dimension for samples

if size(x,1)>1 && size(x,2)==1
    x = x';
end

switch ndims(data)
    case 1, v = interpn(g.vs{1}, data, x);
    case 2, v = interpn(g.vs{1}, g.vs{2}, data, x(:,1),x(:,2));
%     case 3, v = interpn(g.vs{1}, g.vs{2}, [g.vs{3}; 2*pi], cat(3,u,u(:,:,1)),  x(:,1),x(:,2),x(:,3));
    case 3, v = interpn(g.vs{1}, g.vs{2}, g.vs{3}, data,  x(:,1),x(:,2),x(:,3));
    case 4, v = interpn(g.vs{1},g.vs{2},g.vs{3}, g.vs{4}, data, x(:,1),x(:,2),x(:,3),x(:,4));
    case 6
        v = interpn(g.vs{1}, g.vs{2},g.vs{3}, g.vs{4}, g.vs{5}, g.vs{6}, ...
            data, x(:,1), x(:,2), x(:,3), x(:,4), x(:,5), x(:,6));
    otherwise, error(['Cannot evaluate matrices with dimension' num2str(g.dim) '!'])
end

end