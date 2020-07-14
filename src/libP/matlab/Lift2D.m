function [LIFT] = Lift2D(N, faceNodes, r, s)

% function [LIFT] = Lift2D()
% Purpose  : Compute surface to volume lift term for DG formulation
Nfp = N+1;
Np = (N+1)*(N+2)/2;
Nfaces = 3;

Emat = zeros(Np, Nfaces*Nfp);

% face 1
faceR = r(faceNodes(:,1));
V1D = Vandermonde1D(N, faceR); 
massEdge1 = inv(V1D*V1D');
Emat(faceNodes(:,1),1:Nfp) = massEdge1;

% face 2
faceR = r(faceNodes(:,2));
V1D = Vandermonde1D(N, faceR);
massEdge2 = inv(V1D*V1D');
Emat(faceNodes(:,2),Nfp+1:2*Nfp) = massEdge2;

% face 3
faceS = s(faceNodes(:,3));
V1D = Vandermonde1D(N, faceS); 
massEdge3 = inv(V1D*V1D');
Emat(faceNodes(:,3),2*Nfp+1:3*Nfp) = massEdge3;

% inv(mass matrix)*\I_n (L_i,L_j)_{edge_n}
V = Vandermonde2D(N, r, s);
LIFT = V*(V'*Emat);
return
