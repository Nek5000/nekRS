function [xcoords,ycoords] = WarpBlendTransform2D(N, rcoords, scoords, alpha);

% function [xcoords,ycoords] = WarpBlendTransform2D(N, rcoords, scoords);
% Purpose  : Compute Warp & Blend transform (rcoords,scoords)->(xcoords,ycoords) 
%            for nodes in equilateral triangle. Polynomial of order N

if(nargin==3)
  alpopt = [0.0000 0.0000 1.4152 0.1001 0.2751 0.9800 1.0999 ...
	     1.2832 1.3648 1.4773 1.4959 1.5743 1.5770 1.6223 1.6258];

  % Set optimized parameter, alpha, depending on order N
  if (N<16)
    alpha = alpopt(N);
  else
    alpha = 5/3;
  end
end

% save dimension of original coordinates
orig_dim = size(rcoords);

% Create equidistributed nodes on equilateral triangle
L1 =  0.5*(1+scoords(:));
L2 = -0.5*(rcoords(:)+scoords(:));
L3 =  0.5*(1+rcoords(:));

xcoords = -L2+L3; 
ycoords = (-L2-L3+2*L1)/sqrt(3.0);

% Compute blending function at each node for each edge
blend1 = 4*L2.*L3; blend2 = 4*L1.*L3; blend3 = 4*L1.*L2;

% Amount of warp for each node, for each edge
warpf1 = Warpfactor(N,L3-L2); 
warpf2 = Warpfactor(N,L1-L3); 
warpf3 = Warpfactor(N,L2-L1);

% Combine blend & warp
warp1 = blend1.*warpf1.*(1 + (alpha*L1).^2);
warp2 = blend2.*warpf2.*(1 + (alpha*L2).^2);
warp3 = blend3.*warpf3.*(1 + (alpha*L3).^2);

% Accumulate deformations associated with each edge
xcoords = xcoords + 1*warp1 + cos(2*pi/3)*warp2 + cos(4*pi/3)*warp3;
ycoords = ycoords + 0*warp1 + sin(2*pi/3)*warp2 + sin(4*pi/3)*warp3;

% reshape 
xcoords = reshape(xcoords, orig_dim);
ycoords = reshape(ycoords, orig_dim);

return;