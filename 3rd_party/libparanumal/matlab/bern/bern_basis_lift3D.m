function [L0ids,L0vals,ELids,ELvals] = bern_basis_lift3D(N,V,VB,r,s,t)

Np = (N+1)*(N+2)*(N+3)/6;
Nfaces = 4;
Nfp = (N+1)*(N+2)/2;

tol = 1E-3;

[r2D,s2D] = Nodes2D(N); [r2D,s2D] = xytors(r2D,s2D);
VB2D = bern_basis_tri(N,r2D,s2D);

fmask1   = find( abs(t+1) < tol)';
fmask2   = find( abs(s+1) < tol)';
fmask3   = find( abs(r+s+t+1) < tol)';
fmask4   = find( abs(r+1) < tol)';
Fmask  = [fmask1;fmask2;fmask3;fmask4]';

Emat = zeros(Np, Nfaces*Nfp);
% face 1
faceR = r(Fmask(:,1));
faceS = s(Fmask(:,1));
V2D = Vandermonde2D(N, faceR, faceS);
massFace1 = inv(V2D*V2D');
Emat(Fmask(:,1),1:Nfp) = massFace1;
% face 2
faceR = r(Fmask(:,2));
faceT = t(Fmask(:,2));
V2D = Vandermonde2D(N, faceR, faceT);
massFace2 = inv(V2D*V2D');
Emat(Fmask(:,2),Nfp+1:2*Nfp) = massFace2;
% face 3
faceR = r(Fmask(:,3));
faceS = s(Fmask(:,3));
V2D = Vandermonde2D(N, faceR, faceS);
massFace3 = inv(V2D*V2D');
Emat(Fmask(:,3),2*Nfp+1:3*Nfp) = massFace3;
% face 4
faceS = s(Fmask(:,4));
faceT = t(Fmask(:,4));
V2D = Vandermonde2D(N, faceS, faceT);
massFace4 = inv(V2D*V2D');
Emat(Fmask(:,4),3*Nfp+1:4*Nfp) = massFace4;
% inv(mass matrix)*\I_n (L_i,L_j)_{edge_n}
LIFT = V*(V'*Emat);
% convert to BB
LIFT = VB\(LIFT * kron(eye(Nfaces),VB2D));

EL1 = [];
for j = 0:N
    l(j+1) = (-1)^j * nchoosek(N,j)/(1+j);
    Ej = bern_basis_tri(N,r2D,s2D)\bern_basis_tri(N-j,r2D,s2D);
    EL1 = [EL1; l(j+1)*Ej'];
end
EL1(abs(EL1)<tol) = 0;


[r2Dp1 s2Dp1] = Nodes2D(N+1); [r2Dp1 s2Dp1] = xytors(r2Dp1,s2Dp1);
ENp1 = bern_basis_tri(N+1,r2Dp1,s2Dp1)\bern_basis_tri(N,r2Dp1,s2Dp1);

L0 = (N+1)^2/2 * ENp1'*ENp1;
L0(abs(L0)<tol) = 0;

EL = zeros(Np,Nfaces*Nfp);
EL(:,1:Nfp) = EL1;
for f = 2:Nfaces
    ids = (1:Nfp)+(f-1)*Nfp;
    diff = sum(abs( - LIFT(:,ids)),2);
    p = zeros(Np,1);
    for i = 1:Np
        iid = find(sum(abs(repmat(LIFT(i,1:Nfp),Np,1)-LIFT(:,ids)),2)<1e-8);
        p(iid) = i;
    end
    EL(:,ids) = EL1(p,:);
end

EL(abs(EL)<tol) = 0;

% accuracy check
%if norm(LIFT-EL*kron(eye(Nfaces),L0),'fro') > tol
%    keyboard
%end

% L0
L0ids = zeros(Nfp,7);
L0vals = zeros(Nfp,7);
if max(sum(abs(L0)>0,2)) > 7
    keyboard
end
for i = 1:Nfp
    tmp = find(L0(i,:));
    L0vals(i,1:length(tmp)) = L0(i,tmp);
    tmp = tmp-1; % zero indexing
    if length(tmp) < 7
        tmp = [tmp zeros(1,7-length(tmp))];
    end
    L0ids(i,:) = tmp;
end

max_EL_nnz = Nfp + 3;
ELids = zeros(Np,max_EL_nnz);
ELvals = zeros(Np,max_EL_nnz);
if max(sum(abs(EL)>0,2)) ~= max_EL_nnz
    keyboard
end
for i = 1:Np
    tmp = find(EL(i,:));
    ELvals(i,1:length(tmp)) = EL(i,tmp);
    tmp = tmp-1; % zero indexing
    if length(tmp) < max_EL_nnz
        tmp = [tmp zeros(1,max_EL_nnz-length(tmp))];
    end
    ELids(i,:) = tmp;
end


return
