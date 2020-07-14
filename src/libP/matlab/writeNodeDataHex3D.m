
function writeNodeDataHex3D(N, Nc)

Np = (N+1)*(N+1)*(N+1);

r1d = JacobiGL(0,0,N);
V1d = Vandermonde1D(N, r1d);
D1d = Dmatrix1D(N, r1d, V1d);
M1d = inv(V1d')/V1d;
w1d = sum(M1d);
cnt = 1;
r = zeros(Np,1);
s = zeros(Np,1);
t = zeros(Np,1);
for k=1:N+1
  for j=1:N+1
    for i=1:N+1
      r(cnt) = r1d(i); %% r runs fastest
      s(cnt) = r1d(j);
      t(cnt) = r1d(k);
      cnt = cnt+1;
    end
  end
end
Np = (N+1)*(N+1)*(N+1);
r = reshape(r, Np,1);
s = reshape(s, Np,1);
t = reshape(t, Np,1);

Nfp = (N+1)*(N+1);
Nfaces = 6;

% find all the nodes that lie on each edge
NODETOL = 1e-8;
faceNodes = zeros(Nfp, Nfaces);
faceNodes(:,1)   = find( abs(t+1) < NODETOL); 
faceNodes(:,2)   = find( abs(s+1) < NODETOL);
faceNodes(:,3)   = find( abs(r-1) < NODETOL);
faceNodes(:,4)   = find( abs(s-1) < NODETOL);
faceNodes(:,5)   = find( abs(r+1) < NODETOL);
faceNodes(:,6)   = find( abs(t-1) < NODETOL);

V = VandermondeHex3D(N, r, s, t);
[Dr,Ds,Dt] = DmatricesHex3D(N, r, s, t, V);
LIFT = LiftHex3D(N, faceNodes, r, s, t);

fname = sprintf('hexN%02d.dat', N);

fid = fopen(fname, 'w');

writeFloatMatrix(fid, r, 'Nodal r-coordinates');
writeFloatMatrix(fid, s, 'Nodal s-coordinates');
writeFloatMatrix(fid, t, 'Nodal t-coordinates');

if (0) 
writeFloatMatrix(fid, full(Dr), 'Nodal Dr differentiation matrix');
writeFloatMatrix(fid, full(Ds), 'Nodal Ds differentiation matrix');
writeFloatMatrix(fid, full(Dt), 'Nodal Dt differentiation matrix');
writeFloatMatrix(fid, full(LIFT), 'Nodal Lift Matrix');
end

writeIntMatrix(fid, faceNodes'-1, 'Nodal Face nodes');

writeFloatMatrix(fid, r1d, 'Nodal 1D GLL Nodes');
writeFloatMatrix(fid, w1d', 'Nodal 1D GLL Weights');
writeFloatMatrix(fid, D1d, 'Nodal 1D differentiation matrix');

%% compute equispaced nodes on equilateral triangle
[plotR,plotS,plotT] = meshgrid(linspace(-1,1,N+1)); %% hacked
plotR = plotR(:); plotS = plotS(:); plotT = plotT(:);

%% count plot nodes
plotNp = length(plotR);

%% triangulate equilateral element nodes
plotEToV = delaunayFixVolume(plotR,plotS,plotT)-1; 

%% count triangles in plot node triangulation
plotNelements = size(plotEToV,1); 

%% create interpolation matrix from warp & blend to plot nodes
plotInterp = VandermondeHex3D(N, plotR,plotS,plotT)/V; 

writeFloatMatrix(fid, plotR, 'Plotting r-coordinates');
writeFloatMatrix(fid, plotS, 'Plotting s-coordinates');
writeFloatMatrix(fid, plotT, 'Plotting t-coordinates');
writeFloatMatrix(fid, plotInterp, 'Plotting Interpolation Matrix');
writeIntMatrix(fid, plotEToV, 'Plotting triangulation');

%% 1D 
gllS = transpose(D1d)*diag(w1d)*D1d;

NqP = N+3;

%ids 
Nelements = 3;
cnt = 1;
for e=1:Nelements
for n=1:N+1
galnums(n,e) = cnt;
cnt = cnt+1;
end
cnt = cnt-1;
end

A = zeros(cnt,cnt);
for e=1:Nelements
for n=1:N+1
for m=1:N+1
i = galnums(n,e);
j = galnums(m,e);

A(i,j) = A(i,j) + gllS(n,m);
end
end
end

%% WARNING NEED N>1 (otherwise we need a boundary condition)

overlap = 1;
ids = N+1-overlap:2*N+1+overlap;
subA = A(ids,ids);

SP = zeros(NqP,NqP); %% one point overlap
SP(2:NqP-1,2:NqP-1) = gllS;
SP(1,1) = SP(1,1) + gllS(2,2);
SP(2,1) = SP(2,1) + gllS(1,2);
SP(2,2) = SP(2,2) + gllS(1,1);
SP(1,2) = SP(1,2) + gllS(2,1);

SP(NqP,NqP)   = SP(NqP,NqP) + gllS(2,2);
SP(NqP-1,NqP) = SP(NqP-1,NqP) + gllS(1,2);
SP(NqP-1,NqP-1) = SP(NqP-1,NqP-1) + gllS(1,1);
SP(NqP,NqP-1) = SP(NqP,NqP-1) + gllS(2,1);

gllwP = diag([w1d(2),2*w1d(1),w1d(2:end-1),2*w1d(1),w1d(2)]);

[vSP,dSP] = eig(gllwP\SP);

%% invSP = vSP*inv(dSP)*inv(gllwP*vSP) (i.e. the inverse we need)
%% define P = vSP, invP = inv(gllwP*vSP), 
P = vSP;
invP = inv(gllwP*vSP);
diagOp = diag(dSP); % need to divide to get inverse

writeFloatMatrix(fid, invP, 'C0 overlapping patch forward matrix');
writeFloatMatrix(fid, diagOp, 'C0 overlapping patch diagonal scaling');
writeFloatMatrix(fid, P, 'C0 overlapping patch backward matrix');


%%ids 
Nelements = 10;
Nq = N+1;
ADG = zeros(Nq*Nelements,Nq*Nelements);
es = reshape(1:Nelements*Nq, Nq,Nelements);

tau = 2*(N+1)^2;
for e=2:Nelements-1
  n = es(:,e);
  nL = es(1,e);  
  nR = es(Nq,e);
  nP = es(:,e+1);
  nM = es(:,e-1);

  ADG(n,n)  = ADG(n,n)+gllS;

  ADG(n,nL)   = ADG(n,nL)   + 0.5*transpose(D1d(1,:));
  ADG(n,nL-1) = ADG(n,nL-1) - 0.5*transpose(D1d(1,:));

  ADG(n,nR)   = ADG(n,nR)   - 0.5*transpose(D1d(Nq,:));
  ADG(n,nR+1) = ADG(n,nR+1) + 0.5*transpose(D1d(Nq,:));

  ADG(nL,n)   = ADG(nL,n)   + 0.5*(D1d(1,:));
  ADG(nL-1,n) = ADG(nL-1,n) - 0.5*(D1d(1,:));

  ADG(nR,n)   = ADG(nR,n)   - 0.5*(D1d(Nq,:));
  ADG(nR+1,n) = ADG(nR+1,n) + 0.5*(D1d(Nq,:));

  ADG(nL,nL)  = ADG(nL,nL) + 0.5*tau;	    
  ADG(nL,nL-1) = ADG(nL,nL-1) - 0.5*tau;	    

  ADG(nR,nR) = ADG(nR,nR) + 0.5*tau;	    
  ADG(nR,nR+1) = ADG(nR,nR+1) - 0.5*tau;	    

  MDG(n,n) = diag(w1d);
end

ids = 4*Nq:5*Nq+1;
BDG = ADG(ids,ids);
MDG = MDG(ids,ids);

gllwP = diag([w1d(1),w1d,w1d(1)]);

[vSP,dSP] = eig(gllwP\BDG);

%% invSP = vSP*inv(dSP)*inv(gllwP*vSP) (i.e. the inverse we need)
%% define P = vSP, invP = inv(gllwP*vSP), 
P = vSP;
invP = inv(gllwP*vSP);
diagOp = diag(dSP); % need to divide to get inverse

writeFloatMatrix(fid, invP, 'IPDG overlapping patch forward matrix');
writeFloatMatrix(fid, diagOp, 'IPDG overlapping patch diagonal scaling');
writeFloatMatrix(fid, P, 'IPDG overlapping patch backward matrix');


[gr,gw] = JacobiGQ(0,0,N+1);
gI = Vandermonde1D(N, gr)/Vandermonde1D(N, r1d);
gD = GradVandermonde1D(N, gr)/Vandermonde1D(N, r1d);
gNq = length(gr);
gV2 = Vandermonde1D(N+1, gr);
gD2 = Dmatrix1D(N+1, gr, gV2);

writeFloatMatrix(fid, gr, 'Gauss Legendre 1D quadrature nodes');
writeFloatMatrix(fid, gw, 'Gauss Legendre 1D quadrature weights');
writeFloatMatrix(fid, gI, 'GLL to Gauss Legendre interpolation matrix');
writeFloatMatrix(fid, gD, 'GLL to Gauss Legendre differentiation matrix');
writeFloatMatrix(fid, gD2, 'Gauss Legendre to Gauss Legendre differentiation matrix');


%% 1D quadrature
%Nc = ceil(3*N/2);

%[z,w] = JacobiGQ(0,0,Nc);
[z] = JacobiGL(0,0,Nc);
w = sum(inv(Vandermonde1D(Nc, z)*Vandermonde1D(Nc,z)'))'


[N, Nc, length(w)]
cInterp = Vandermonde1D(N, z)/Vandermonde1D(N, r1d);
cubProject = (cInterp)';
cubD = Dmatrix1D(Nc, z, Vandermonde1D(Nc,z));
cubDT = (Dmatrix1D(Nc, z, Vandermonde1D(Nc,z)))';
Nqc = length(z)

writeFloatMatrix(fid, z, 'Quadrature r-coordinates');
writeFloatMatrix(fid, w, 'Quadrature weights');

writeFloatMatrix(fid, cInterp, 'Quadrature Interpolation Matrix');
writeFloatMatrix(fid, cubDT, 'Quadrature Weak D Differentiation Matrix');
writeFloatMatrix(fid, cubD, 'Quadrature Differentiation Matrix');
writeFloatMatrix(fid, cubProject, 'Quadrature Projection Matrix');


%% volume cubature
%{

[z,w] = JacobiGQ(0,0,ceil(3*N/2));
Nz = length(z);
sk = 1;
cubr = zeros(Nz*Nz*Nz,1);
cubs = zeros(Nz*Nz*Nz,1);
cubt = zeros(Nz*Nz*Nz,1);
cubw = zeros(Nz*Nz*Nz,1);
for k=1:Nz
for j=1:Nz
for i=1:Nz
cubr(sk) = z(i);
cubs(sk) = z(j);
cubt(sk) = z(k);
cubw(sk) = w(i)*w(j)*w(k);
sk = sk+1;
end
end
end

cInterp = VandermondeHex3D(N, cubr, cubs, cubt)/V;
Ncub = length(cubr);

cV = VandermondeHex3D(N, cubr, cubs, cubt);
cV'*diag(cubw)*cV;

[cVr,cVs,cVt] = GradVandermondeHex3D(N, cubr, cubs, cubt);
cubDrT = V*transpose(cVr)*diag(cubw);
cubDsT = V*transpose(cVs)*diag(cubw);
cubDtT = V*transpose(cVt)*diag(cubw);
cubProject = V*cV'*diag(cubw); %% relies on (transpose(cV)*diag(cubw)*cV being the identity)

writeFloatMatrix(fid, cubr, 'Cubature r-coordinates');
writeFloatMatrix(fid, cubs, 'Cubature s-coordinates');
writeFloatMatrix(fid, cubt, 'Cubature t-coordinates');
writeFloatMatrix(fid, cubw, 'Cubature weights');

writeFloatMatrix(fid, cInterp, 'Cubature Interpolation Matrix');
writeFloatMatrix(fid, cubDrT, 'Cubature Weak Dr Differentiation Matrix');
writeFloatMatrix(fid, cubDsT, 'Cubature Weak Ds Differentiation Matrix');
writeFloatMatrix(fid, cubDtT, 'Cubature Weak Dt Differentiation Matrix');
writeFloatMatrix(fid, cubProject, 'Cubature Projection Matrix');

end
%}

%degree raising interpolation
rP1 = JacobiGL(0,0,N+1);
VP1 = Vandermonde1D(N, rP1);

IP1 = VP1/V1d;
NpP1 = length(rP1);

%degree lowering interpolation
if(N>1)
  rM1 = JacobiGL(0,0,N-1);
else %hard code degree 0
  rM1 = 0;
end

VM1 = Vandermonde1D(N, rM1);
IM1 = VM1/V1d;
NpM1 = length(rM1);

writeFloatMatrix(fid, IP1, '1D degree raise matrix');
writeFloatMatrix(fid, IM1, '1D degree lower matrix');


%%SEMFEM data
%manually build quad grid for FEM problem
FEMEToV = zeros(N*N*N,8);
cnt =1;
for k=0:N-1
  for j=0:N-1
    for i=0:N-1
      id = i+(j)*(N+1)+k*(N+1)*(N+1);
      FEMEToV(cnt,:) = [id, id+1, id+1+(N+1), id+(N+1), id+(N+1)*(N+1), id+1+(N+1)*(N+1), id+1+(N+1)+(N+1)*(N+1), id+(N+1)+(N+1)*(N+1)];
      cnt = cnt+1;
    end
  end
end

writeIntMatrix(fid, FEMEToV, 'SEMFEM reference mesh');

fclose(fid);

