function [r,s] = xytors(x,y)

% function [r,s] = xytors(x, y)
% Purpose : Transfer from (x,y) in equilateral triangle
%           to (r,s) coordinates in standard triangle

L1 = (sqrt(3.0)*y+1.0)/3.0;
L2 = (-3.0*x - sqrt(3.0)*y + 2.0)/6.0;
L3 = ( 3.0*x - sqrt(3.0)*y + 2.0)/6.0;

r = -L2 + L3 - L1; s = -L2 - L3 + L1;
return;
