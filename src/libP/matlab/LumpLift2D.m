function [LIFT] = LumpLift2D()

% function [LIFT] = Lift2D()
% Purpose  : Compute surface to volume lift term for DG formulation

Globals2D;
Emat = zeros(Np, Nfaces*Nfp);

% face 1
faceR = r(Fmask(:,1));
V1D = Vandermonde1D(N, faceR); 
massEdge1 = diag(sum(inv(V1D*V1D')))
Emat(Fmask(:,1),1:Nfp) = massEdge1;

% face 2
faceR = r(Fmask(:,2));
V1D = Vandermonde1D(N, faceR);
massEdge2 = diag(sum(inv(V1D*V1D')))
Emat(Fmask(:,2),Nfp+1:2*Nfp) = massEdge2;

% face 3
faceS = s(Fmask(:,3));
V1D = Vandermonde1D(N, faceS); 
massEdge3 = diag(sum(inv(V1D*V1D')))
Emat(Fmask(:,3),2*Nfp+1:3*Nfp) = massEdge3;

% inv(mass matrix)*\I_n (L_i,L_j)_{edge_n}
LIFT = V*(V'*Emat);
return
