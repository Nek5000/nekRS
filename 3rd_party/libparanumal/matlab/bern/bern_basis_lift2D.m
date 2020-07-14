function [L0vals ELids ELvals] = bern_basis_lift2D(N,V,VB,r,s)

Np = (N+1)*(N+2)/2;
Nfaces = 3;
Nfp = N+1;

tol = 1E-3;

fmask1   = find( abs(s+1) < tol)';
fmask2   = find( abs(r+s) < tol)';
fmask3   = find( abs(r+1) < tol)';
Fmask  = [fmask1;fmask2;fmask3]';

[r1Dq w1Dq] = JacobiGQ(0,0,N);

Emat = zeros(Np, Nfaces*Nfp);
VBmat = zeros(Nfaces*Nfp,Nfaces*Nfp);

% face 1
faceR = r(Fmask(:,1));
V1D = Vandermonde1D(N, faceR);
massEdge1 = inv(V1D*V1D');
Emat(Fmask(:,1),1:Nfp) = massEdge1;
VB1D = bern_basis_1D(N,faceR);
VBmat(1:Nfp,1:Nfp) = VB1D;

% face 2
faceR = r(Fmask(:,2));
V1D = Vandermonde1D(N, faceR);
massEdge2 = inv(V1D*V1D');
Emat(Fmask(:,2),Nfp+1:2*Nfp) = massEdge2;
VB1D = bern_basis_1D(N,faceR);
VBmat(Nfp+1:2*Nfp,Nfp+1:2*Nfp) = VB1D;

% face 3
faceS = s(Fmask(:,3));
V1D = Vandermonde1D(N, faceS);
massEdge3 = inv(V1D*V1D');
Emat(Fmask(:,3),2*Nfp+1:3*Nfp) = massEdge3;
VB1D = bern_basis_1D(N,faceS);
VBmat(2*Nfp+1:3*Nfp,2*Nfp+1:3*Nfp) = VB1D;

EL1 = [];
for j = 0:N
    l(j+1) = (-1)^j * nchoosek(N,j)/(1+j);
    Ej = bern_basis_1D(N,r1Dq)\bern_basis_1D(N-j,r1Dq);
    EL1 = [EL1; l(j+1)*Ej'];
end
EL1(abs(EL1)<tol) = 0;

ENp1 = bern_basis_1D(N+1,JacobiGL(0,0,N+1))\bern_basis_1D(N,JacobiGL(0,0,N+1));

L0 = (N+1)^2/2 * ENp1'*ENp1;

LIFT = V*(V'*Emat);
% convert to BB
LIFT = VB\(LIFT * VBmat);

EL = zeros(Np,Nfaces*Nfp);
EL(:,1:Nfp) = EL1;

%fill the remaining faces of EL by finding permutations of face1
for f = 2:Nfaces
    ids = (1:Nfp)+(f-1)*Nfp;
    p = zeros(Np,1);
    for i = 1:Np
        iid = find(sum(abs(repmat(LIFT(i,1:Nfp),Np,1)-LIFT(:,ids)),2)<tol);
        p(iid) = i;
    end
    EL(:,ids) = EL1(p,:);
end

EL(abs(EL)<tol) = 0;

% L0 vals (L0 is tridiagonal)
L0vals = zeros(Nfp,3);
for i = 1:Nfp
    if (i==1)
        L0vals(i,2:3) = [L0(i,i) L0(i,i+1)]; %tridiagonal fix
    elseif (i==Nfp)
        L0vals(i,1:2) = [L0(i,i-1) L0(i,i)];
    else
        L0vals(i,:) = [L0(i,i-1) L0(i,i) L0(i,i+1)];
    end
end

max_EL_nnz = Nfp + 2;
ELids = zeros(Np,max_EL_nnz);
ELvals = zeros(Np,max_EL_nnz);

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
