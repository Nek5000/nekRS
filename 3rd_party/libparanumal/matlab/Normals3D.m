function [nx, ny, nz, sJ] = Normals3D()

% function [nx, ny, nz, sJ] = Normals3D()
% Purpose : Compute outward pointing normals at elements faces as well as surface Jacobians

Globals3D;

[rx,sx,tx,ry,sy,ty,rz,sz,tz,J] = GeometricFactors3D(x,y,z,Dr,Ds,Dt);

% interpolate geometric factors to face nodes
frx = rx(Fmask(:), :); fsx = sx(Fmask(:), :); ftx = tx(Fmask(:), :);
fry = ry(Fmask(:), :); fsy = sy(Fmask(:), :); fty = ty(Fmask(:), :);
frz = rz(Fmask(:), :); fsz = sz(Fmask(:), :); ftz = tz(Fmask(:), :);

% build normals
nx = zeros(4*Nfp, K); ny = zeros(4*Nfp, K); nz = zeros(4*Nfp, K);
fid1 = (1:Nfp)'; fid2 = (Nfp+1:2*Nfp)'; 
fid3 = (2*Nfp+1:3*Nfp)'; fid4 = (3*Nfp+1:4*Nfp)';

% face 1
nx(fid1, :) = -ftx(fid1, :); ny(fid1, :) = -fty(fid1, :); nz(fid1, :) = -ftz(fid1, :);

% face 2
nx(fid2, :) = -fsx(fid2, :); ny(fid2, :) = -fsy(fid2, :); nz(fid2, :) = -fsz(fid2, :);

% face 3
nx(fid3, :) = frx(fid3, :)+fsx(fid3, :)+ftx(fid3, :);
ny(fid3, :) = fry(fid3, :)+fsy(fid3, :)+fty(fid3, :);
nz(fid3, :) = frz(fid3, :)+fsz(fid3, :)+ftz(fid3, :);

% face 4
nx(fid4, :) = -frx(fid4, :); ny(fid4, :) = -fry(fid4, :); nz(fid4, :) = -frz(fid4, :);

% normalise
sJ = sqrt(nx.*nx+ny.*ny+nz.*nz);
nx = nx./sJ; ny = ny./sJ; nz = nz./sJ;
sJ = sJ.*J(Fmask(:), :);
return;
