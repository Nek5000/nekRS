function [nx, ny, sJ] = Normals2D()

% function [nx, ny, sJ] = Normals2D()
% Purpose : Compute outward pointing normals at
%	    elements faces as well as surface Jacobians

Globals2D;
xr = Dr*x; yr = Dr*y; xs = Ds*x; ys = Ds*y; J = xr.*ys-xs.*yr;

% interpolate geometric factors to face nodes
fxr = xr(Fmask, :); fxs = xs(Fmask, :);
fyr = yr(Fmask, :); fys = ys(Fmask, :);

% build normals
nx = zeros(3*Nfp, K); ny = zeros(3*Nfp, K);
fid1 = (1:Nfp)'; fid2 = (Nfp+1:2*Nfp)'; fid3 = (2*Nfp+1:3*Nfp)';

% face 1
nx(fid1, :) =  fyr(fid1, :); ny(fid1, :) = -fxr(fid1, :);

% face 2
nx(fid2, :) =  fys(fid2, :)-fyr(fid2, :);
ny(fid2, :) = -fxs(fid2, :)+fxr(fid2, :);

% face 3
nx(fid3, :) = -fys(fid3, :); ny(fid3, :) =  fxs(fid3, :);

% normalise
sJ = sqrt(nx.*nx+ny.*ny); nx = nx./sJ; ny = ny./sJ;
return;
