function [LIFT] = LiftHex3D(N, faceNodes, r, s, t)

% function [LIFT] = LiftQuad3D()
% Purpose  : Compute surface to volume lift term for DG formulation

Nfp = (N+1)*(N+1);
Nfaces = 6;
Np = length(r);

Emat = zeros(Np, Nfaces*Nfp);

% face 1
faceR = r(faceNodes(:,1));
faceS = s(faceNodes(:,1));
V2D = VandermondeQuad2D(N, faceR, faceS); 
massEdge1 = inv(V2D*V2D');
Emat(faceNodes(:,1),1:Nfp) = massEdge1;

% face 2
faceR = r(faceNodes(:,2));
faceS = t(faceNodes(:,2));
V2D = VandermondeQuad2D(N, faceR, faceS); 
massEdge2 = inv(V2D*V2D');
Emat(faceNodes(:,2),Nfp+1:2*Nfp) = massEdge2;

% face 3
faceR = s(faceNodes(:,3));
faceS = t(faceNodes(:,3));
V2D = VandermondeQuad2D(N, faceR, faceS); 
massEdge3 = inv(V2D*V2D');
Emat(faceNodes(:,3),2*Nfp+1:3*Nfp) = massEdge3;

% face 4
faceR = r(faceNodes(:,4));
faceS = t(faceNodes(:,4));
V2D = VandermondeQuad2D(N, faceR, faceS); 
massEdge4 = inv(V2D*V2D');
Emat(faceNodes(:,4),3*Nfp+1:4*Nfp) = massEdge4;

% face 5
faceR = s(faceNodes(:,5));
faceS = t(faceNodes(:,5));
V2D = VandermondeQuad2D(N, faceR, faceS); 
massEdge5 = inv(V2D*V2D');
Emat(faceNodes(:,5),4*Nfp+1:5*Nfp) = massEdge5;

% face 6
faceR = r(faceNodes(:,6));
faceS = s(faceNodes(:,6));
V2D = VandermondeQuad2D(N, faceR, faceS); 
massEdge6 = inv(V2D*V2D');
Emat(faceNodes(:,6),5*Nfp+1:6*Nfp) = massEdge6;

% inv(mass matrix)*\I_n (L_i,L_j)_{edge_n}
V = VandermondeHex3D(N, r, s, t);
LIFT = V*(V'*Emat);

ids = find(abs(LIFT)<1e-13); LIFT(ids) = 0;
LIFT = sparse(LIFT);

return
