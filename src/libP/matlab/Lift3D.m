function [LIFT] = Lift3D(N, faceNodes, r, s, t)

% function [LIFT] = Lif3D()
% Purpose  : Compute surface to volume lift term for DG formulation
  Np = length(r);
Nfp = (N+1)*(N+2)/2;
Nfaces = 4;

Emat = zeros(Np, Nfaces*Nfp);

% face 1
faceR = r(faceNodes(:,1));
faceS = s(faceNodes(:,1));
V2D = Vandermonde2D(N, faceR, faceS); 
massEdge1 = inv(V2D*V2D');
Emat(faceNodes(:,1),1:Nfp) = massEdge1;

% face 2
faceR = r(faceNodes(:,2));
faceT = t(faceNodes(:,2));
V2D = Vandermonde2D(N, faceR, faceT);
massEdge2 = inv(V2D*V2D');
Emat(faceNodes(:,2),Nfp+1:2*Nfp) = massEdge2;

% face 3
faceS = s(faceNodes(:,3));
faceT = t(faceNodes(:,3));
V2D = Vandermonde2D(N, faceS, faceT); 
massEdge3 = inv(V2D*V2D');
Emat(faceNodes(:,3),2*Nfp+1:3*Nfp) = massEdge3;

% face 4
faceS = s(faceNodes(:,4));
faceT = t(faceNodes(:,4));
V2D = Vandermonde2D(N, faceS, faceT); 
massEdge4 = inv(V2D*V2D');
Emat(faceNodes(:,4),3*Nfp+1:4*Nfp) = massEdge4;

% inv(mass matrix)*\I_n (L_i,L_j)_{edge_n}
V = Vandermonde3D(N, r, s, t);
LIFT = V*(V'*Emat);
return
