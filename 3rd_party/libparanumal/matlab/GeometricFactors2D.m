function [rx,sx,ry,sy,J] = GeometricFactors2D(x,y,Dr,Ds)

% function [rx,sx,ry,sy,J] = GeometricFactors2D(x,y,Dr,Ds)
% Purpose  : Compute the metric elements for the local
%            mappings of the elements

% Calculate geometric factors
xr = Dr*x; xs = Ds*x; yr = Dr*y; ys = Ds*y; J = -xs.*yr + xr.*ys;
rx = ys./J; sx =-yr./J; ry =-xs./J; sy = xr./J;
return;
