function [LIFT] = LiftQuad2D(N, faceNodes, r, s)

% function [LIFT] = LiftQuad2D()
% Purpose  : Compute surface to volume lift term for DG formulation

Nfp = N+1;
Nfaces = 4;
Np = length(r);

Emat = zeros(Np, Nfaces*Nfp);

% face 1
faceR = r(faceNodes(:,1));
V1D = Vandermonde1D(N, faceR); 
massEdge1 = inv(V1D*V1D');
Emat(faceNodes(:,1),1:Nfp) = massEdge1;

% face 2
faceR = s(faceNodes(:,2));
V1D = Vandermonde1D(N, faceR);
massEdge2 = inv(V1D*V1D');
Emat(faceNodes(:,2),Nfp+1:2*Nfp) = massEdge2;

% face 3
faceR = r(faceNodes(:,3));
V1D = Vandermonde1D(N, faceR); 
massEdge3 = inv(V1D*V1D');
Emat(faceNodes(:,3),2*Nfp+1:3*Nfp) = massEdge3;

% face 4
faceS = s(faceNodes(:,4));
V1D = Vandermonde1D(N, faceS); 
massEdge4 = inv(V1D*V1D');
Emat(faceNodes(:,4),3*Nfp+1:4*Nfp) = massEdge4;

% inv(mass matrix)*\I_n (L_i,L_j)_{edge_n}
V = VandermondeQuad2D(N, r, s);
LIFT = V*(V'*Emat);

ids = find(abs(LIFT)<1e-13); LIFT(ids) = 0;
LIFT = sparse(LIFT);

return
