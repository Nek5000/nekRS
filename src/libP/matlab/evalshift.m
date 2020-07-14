function [dx, dy] = evalshift(p, pval, L1, L2, L3)  

% function [dx, dy] = evalshift(p, pval, L1, L2, L3)  
% Purpose: compute two-dimensional Warp & Blend transform

% 1) compute Gauss-Lobatto-Legendre node distribution
gaussX = -JacobiGL(0,0,p);
 
% 2) compute blending function at each node for each edge
blend1 = L2.*L3; blend2 = L1.*L3; blend3 = L1.*L2;

% 3) amount of warp for each node, for each edge
warpfactor1 = 4*evalwarp(p, gaussX, L3-L2); 
warpfactor2 = 4*evalwarp(p, gaussX, L1-L3); 
warpfactor3 = 4*evalwarp(p, gaussX, L2-L1); 

% 4) combine blend & warp
warp1 = blend1.*warpfactor1.*(1 + (pval*L1).^2);
warp2 = blend2.*warpfactor2.*(1 + (pval*L2).^2);
warp3 = blend3.*warpfactor3.*(1 + (pval*L3).^2);

% 5) evaluate shift in equilateral triangle
dx = 1*warp1 + cos(2*pi/3)*warp2 + cos(4*pi/3)*warp3;
dy = 0*warp1 + sin(2*pi/3)*warp2 + sin(4*pi/3)*warp3;
return;
