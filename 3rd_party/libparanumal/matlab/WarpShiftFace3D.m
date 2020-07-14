function [warpx, warpy] = WarpShiftFace3D(p,pval, pval2, L1,L2,L3,L4)

% function [warpx, warpy] = WarpShiftFace3D(p,pval, pval2, L1,L2,L3,L4)     
% Purpose: compute warp factor used in creating 3D Warp & Blend nodes

[dtan1,dtan2] = evalshift(p, pval, L2, L3, L4);
warpx = dtan1; warpy = dtan2;
return;
