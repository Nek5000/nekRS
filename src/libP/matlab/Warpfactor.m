function warp = Warpfactor(N, rout)

% function warp = Warpfactor(N, rout)
% Purpose  : Compute scaled warp function at order N 
%            based on rout interpolation nodes

% Compute LGL and equidistant node distribution
LGLr = JacobiGL(0,0,N); req  = linspace(-1,1,N+1)';

% Compute V based on req
Veq = Vandermonde1D(N,req);

% Evaluate Lagrange polynomial at rout
Nr = length(rout); Pmat = zeros(N+1,Nr);
for i=1:N+1
  Pmat(i,:) = JacobiP(rout, 0, 0, i-1)';
end;
Lmat = Veq'\Pmat;

% Compute warp factor
warp = Lmat'*(LGLr - req);

% Scale factor
zerof = (abs(rout)<1.0-1.0e-10); sf = 1.0 - (zerof.*rout).^2;
warp = warp./sf + warp.*(zerof-1);
return;
