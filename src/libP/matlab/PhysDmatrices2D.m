function [Dx,Dy] = PhysDmatrices2D(x1, y1, interp)

% function [Dr,Ds] = PhysDmatrices2D(x1, y1, interp)
% Purpose : Initialize the (x,y) differentiation matrices
%	    on the simplex, evaluated at (x1,y1) at order N

Globals2D;
IDr = interp*Dr; IDs = interp*Ds;

[rx1,sx1,ry1,sy1,J1] = GeometricFactors2D(x1, y1, IDr, IDs);

n = size(interp, 1);
Dx = spdiags(rx1, 0, n, n)*IDr + spdiags(sx1, 0, n, n)*IDs;
Dy = spdiags(ry1, 0, n, n)*IDr + spdiags(sy1, 0, n, n)*IDs;

return;
