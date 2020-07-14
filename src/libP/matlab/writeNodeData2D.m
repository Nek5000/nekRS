
function writeNodeData2D(inN)

Globals2D;

N = inN;
Nfp = N+1;
Nfaces = 3;

%% Nodal data
[r,s] = Nodes2D(N);
[r,s] = xytors(r,s);
Np = length(r);

% find all the nodes that lie on each edge
NODETOL = 1e-8;
faceNodes1   = find( abs(s+1) < NODETOL)';
faceNodes2   = find( abs(r+s) < NODETOL)';
faceNodes3   = find( abs(r+1) < NODETOL)';
FaceNodes  = [faceNodes1;faceNodes2;faceNodes3]';

V = Vandermonde2D(N, r, s);
MM = inv(transpose(V))/V;
[Dr,Ds] = Dmatrices2D(N, r, s, V);
LIFT = Lift2D(N, FaceNodes, r, s);

fname = sprintf('triangleN%02d.dat', N);
fid = fopen(fname, 'w');

writeFloatMatrix(fid, r, 'Nodal r-coordinates');
writeFloatMatrix(fid, s, 'Nodal s-coordinates');
writeFloatMatrix(fid, Dr, 'Nodal Dr differentiation matrix');
writeFloatMatrix(fid, Ds, 'Nodal Ds differentiation matrix');
writeFloatMatrix(fid, MM, 'Nodal Mass Matrix');

writeIntMatrix(fid, FaceNodes'-1, 'Nodal Face nodes');
writeFloatMatrix(fid, LIFT, 'Nodal Lift Matrix');

%% Plotting data
%compute equispaced nodes on equilateral triangle
[plotR,plotS] = EquiNodes2D(N+4);
plotNp = length(plotR);
plotEToV = delaunay(plotR,plotS)-1;
plotNelements = size(plotEToV,1);
[plotR,plotS] = xytors(plotR,plotS);

%check triangulation
before = plotNelements;
sk = 0;
for e=1:plotNelements
  v1 = plotEToV(e,1)+1;
  v2 = plotEToV(e,2)+1;
  v3 = plotEToV(e,3)+1;

  x1 = plotR(v1);
  x2 = plotR(v2);
  x3 = plotR(v3);

  y1 = plotS(v1);
  y2 = plotS(v2);
  y3 = plotS(v3);

  plotA = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1);
  if(abs(plotA)>1e-5) sk = sk+1; plotEToV(sk,:) = [v1-1,v2-1,v3-1]; end
end
plotNelements = sk;
plotEToV = plotEToV(1:sk,:);
after = plotNelements;

plotInterp = Vandermonde2D(N, plotR,plotS)/V;

writeFloatMatrix(fid, plotR, 'Plotting r-coordinates');
writeFloatMatrix(fid, plotS, 'Plotting s-coordinates');
writeFloatMatrix(fid, plotInterp, 'Plotting Interpolation Matrix');
writeIntMatrix(fid, plotEToV, 'Plotting triangulation');

%% Cubature data
%volume cubature
[cubr,cubs,cubw] = Cubature2D(3*N);
cInterp = Vandermonde2D(N, cubr, cubs)/V;
Ncub = length(cubr);

cV = Vandermonde2D(N, cubr, cubs);
cV'*diag(cubw)*cV;

[cVr,cVs] = GradVandermonde2D(N, cubr, cubs);
cubDrT = V*transpose(cVr)*diag(cubw);
cubDsT = V*transpose(cVs)*diag(cubw);
cubProject = V*cV'*diag(cubw); %% relies on (transpose(cV)*diag(cubw)*cV being the identity)

writeFloatMatrix(fid, cubr, 'Cubature r-coordinates');
writeFloatMatrix(fid, cubs, 'Cubature s-coordinates');
writeFloatMatrix(fid, cubw, 'Cubature weights');

writeFloatMatrix(fid, cInterp, 'Cubature Interpolation Matrix');
writeFloatMatrix(fid, cubDrT, 'Cubature Weak Dr Differentiation Matrix');
writeFloatMatrix(fid, cubDsT, 'Cubature Weak Ds Differentiation Matrix');
writeFloatMatrix(fid, cubProject, 'Cubature Projection Matrix');

%surface cubature
[z,w] = JacobiGQ(0,0,ceil(3*N/2));
Nfi = length(z);

ir = [z,-z,-ones(Nfi,1)];
is = [-ones(Nfi,1), z, -z];
iw = [w,w,w];

sV = Vandermonde2D(N, ir(:), is(:));
sInterp = sV/V;
interp = [sInterp(1:Nfi,FaceNodes(:,1));sInterp(Nfi+1:2*Nfi,FaceNodes(:,2));sInterp(2*Nfi+1:3*Nfi,FaceNodes(:,3))];

%integration node lift matrix
iLIFT = V*V'*sInterp'*diag(iw(:));

writeFloatMatrix(fid, interp, 'Cubature Surface Interpolation Matrix');
writeFloatMatrix(fid, iLIFT, 'Cubature Surface Lift Matrix');


%% Berstein-Bezier operators 
addpath('./bern')
tol=1e-5;
r1D = JacobiGL(0,0,N);
[r1Dq,w1Dq] = JacobiGQ(0,0,N);
V1Dq = bern_basis_1D(N,r1Dq);
VB1D = bern_basis_1D(N,r1D);

[VB,Vr,Vs,V0,V1,V2] = bern_basis_tri(N,r,s);
[D0ids,D0vals,D1ids,D1vals,D2ids,D2vals] = bern_basis_diff2D(N,V,VB,V0,V1,V2);
[L0vals,ELids,ELvals] = bern_basis_lift2D(N,V,VB,r,s);

invVB = inv(VB);
keyboard
%write out the BB operators
writeFloatMatrix(fid, VB, 'Bernstein-Bezier Vandermonde Matrix');
writeFloatMatrix(fid, invVB, 'Bernstein-Bezier Inverse Vandermonde Matrix');
writeIntMatrix(fid, D0ids, 'Bernstein-Bezier sparse D1 differentiation ids');
writeIntMatrix(fid, D1ids, 'Bernstein-Bezier sparse D2 differentiation ids');
writeIntMatrix(fid, D2ids, 'Bernstein-Bezier sparse D3 differentiation ids');
writeFloatMatrix(fid, D0vals, 'Bernstein-Bezier sparse D differentiation values');

writeFloatMatrix(fid, L0vals, 'Bernstein-Bezier L0 Matrix values');
writeIntMatrix(fid, ELids, 'Bernstein-Bezier EL lift ids');
writeFloatMatrix(fid, ELvals, 'Bernstein-Bezier EL lift values');

%BB cubature
[rq sq wq] = Cubature2D(3*N);
VBq = bern_basis_tri(N,rq,sq);
M = VBq'*diag(wq)*VBq;
PBq = M\(VBq'*diag(wq));

writeFloatMatrix(fid, VBq, 'Cubature Bernstein-Bezier Interpolation Matrix');
writeFloatMatrix(fid, PBq, 'Cubature Bernstein-Bezier Projection Matrix');

%1D degree raise/lower operators
[r1Dq w1Dq] = JacobiGQ(0,0,N+1);
BBRaise = bern_basis_1D(N+1,r1Dq)\bern_basis_1D(N,r1Dq);
BBRaise(abs(BBRaise)<tol) = 0;

BBRaiseIds  = zeros(N+2,2);
BBRaiseVals = zeros(N+2,2);
for i = 1:N+2
  tmp = find(BBRaise(i,:));
  BBRaiseVals(i,1:length(tmp)) = BBRaise(i,tmp);
  tmp = tmp-1; % zero indexing
  if length(tmp) < 2
    tmp = [tmp zeros(1,2-length(tmp))];
  end
  BBRaiseIds(i,:) = tmp;
end

[r1D] = JacobiGQ(0,0,N);
VB1D = bern_basis_1D(N,r1D);
V1D = Vandermonde1D(N, r1D);

[r1Dm1] = JacobiGQ(0,0,N-1);
VB1Dm1 = bern_basis_1D(N-1,r1Dm1);
V1Dm1 = Vandermonde1D(N-1, r1Dm1);

BBLower = V1D\VB1D;
BBLower = VB1Dm1\V1Dm1*BBLower(1:N,1:N+1);

writeIntMatrix(fid, BBRaiseIds, 'Bernstein-Bezier sparse 1D degree raise ids');
writeFloatMatrix(fid, BBRaiseVals, 'Bernstein-Bezier sparse 1D degree raise values');
writeFloatMatrix(fid, BBLower, 'Bernstein-Bezier sparse 1D degree lower matrix');


%% IPDG elliptic patch problem
K = 4;

%VX = [-1,+1,-1,+1,+1,-3];
%VY = [-1,-1,+1,-3,+1,+1];

VX = [-1,1,0,0,2,-2];
VY = [0,0,sqrt(3),-sqrt(3),sqrt(3),sqrt(3)];

EToV = [1,2,3;
	2,1,4;
	3,2,5;
	1,3,6];

	      bc = 3; %Dirichlet;
BCType = [0,0,0;
	  0,bc,bc;
	  0,bc,bc;
	  0,bc,bc]; %3 is a special flag for unconnected internal edge

StartUp2D;

% choose order to integrate exactly
Nint = ceil(2*N/2);

% build cubature nodes for all elements
CubatureOrder = 2*(Nint+1); 
cub = CubatureVolumeMesh2D(CubatureOrder);
  
% build Gauss node data for all element faces
NGauss = (Nint+1); 
gauss = GaussFaceMesh2D(NGauss);
  
% build weak Poisson operator matrices
[A, M] = CurvedPoissonIPDG2D();

full(A);

%% hack since we know W&B face 1 nodes are first
vmapP = reshape(vmapP, Nfp*Nfaces, K);
idsP = vmapP(:,1);
subind = [(1:Np)';idsP];

subA = full(A(subind,subind));
subM = full(M(subind,subind));

condSubA = cond(subA);

[B,d] = eig(subA, subM);

%% A = S + lambda*M 
%%   = M*(M\S + lambda*I) 
%%   ~ J*Mhat*(Mhat\Shat/hscale2 + lambda*I) 
%%   ~ J*Mhat*Bhat*(d/scale + lambda*I)/Bhat
%% inv(A) ~ Bhat*inv(J*(d/scale+lambda*I))*inv(Mhat*Bhat)

forwardMatrix = inv(subM*B);
diagOp = diag(d);
backwardMatrix = B;

NpP = size(subA,1);

writeFloatMatrix(fid, forwardMatrix, 'IPDG overlapping patch forward matrix');
writeFloatMatrix(fid, diagOp, 'IPDG overlapping patch diagonal scaling');
writeFloatMatrix(fid, backwardMatrix, 'IPDG overlapping patch backward matrix');


[rG,sG,shiftIds] = GroupNodes2D(N);
A = full(A);
invA = inv(A);

writeIntMatrix(fid, shiftIds-1, 'Nodal rotation permutations');

writeFloatMatrix(fid, invA, 'IPDG full reference patch inverse matrix');

%degree raising interpolation
[rP1,sP1] = Nodes2D(N+1);
[rP1,sP1] = xytors(rP1,sP1);

VP1 = Vandermonde2D(N, rP1, sP1);
IP1 = VP1/V;
NpP1 = length(rP1);

%degree lowering interpolation
if(N>1)
  [rM1,sM1] = Nodes2D(N-1);
  [rM1,sM1] = xytors(rM1,sM1);
else %hard code degree 0
  rM1 = -1/3;
  sM1 = -1/3;
end

VM1 = Vandermonde2D(N, rM1, sM1);
IM1 = VM1/V;
NpM1 = length(rM1);

writeFloatMatrix(fid, IP1, 'Nodal degree raise matrix');
writeFloatMatrix(fid, IM1, 'Nodal degree lower matrix');

%% SEMFEM
if (N<13)
  addpath('./newNodes')
  [req,seq] = NewEquiNodes2D(N+1,'EI');
  FEMEToV = FemEToV2D(N+1,req,seq,'EI')-1;
  [rFEM,sFEM] = NewNodes2D(N,'EIKappaNp1');
  [rFEM,sFEM] = xytors(rFEM,sFEM);

  FEMEToV
  
  req
  seq
  
  rFEM
  sFEM
  
  triplot(FEMEToV+1,req,seq)

  NpFEM = length(rFEM);
  NelFEM = size(FEMEToV,1);

  IQN = Vandermonde2D(N, rFEM, sFEM)/V;
  invIQN = (transpose(IQN)*IQN)\(transpose(IQN));

  writeFloatMatrix(fid, rFEM, 'SEMFEM r-coordinates');
  writeFloatMatrix(fid, sFEM, 'SEMFEM s-coordinates');

  writeIntMatrix(fid, FEMEToV, 'SEMFEM reference mesh');  
  writeFloatMatrix(fid, invIQN', 'SEMFEM interpolation matrix');
end

%% Sparse Basis
addpath('./sparseBasis')
[cV,cMM,cSrr,cSrs,cSss,stackedNz] = GenModalOps(N);

faceModes1   = find( sum(abs(cV(faceNodes1,:)),1) > NODETOL);
faceModes2   = find( sum(abs(cV(faceNodes2,:)),1) > NODETOL);
faceModes3   = find( sum(abs(cV(faceNodes3,:)),1) > NODETOL);
FaceModes  = [faceModes1;faceModes2;faceModes3]';

for n=1:Np
  for m=1:size(stackedNz,2)
    id = stackedNz(n,m);
    if (id>0)
      sparseSrr(n,m) = cSrr(id,n);
      sparseSrs(n,m) = cSrs(id,n);
      sparseSss(n,m) = cSss(id,n);
     end
   end
end

writeFloatMatrix(fid, cV, 'Sparse basis Vandermonde');
writeFloatMatrix(fid, cMM, 'Sparse basis mass matrix');

writeIntMatrix(fid, FaceModes'-1, 'Sparse basis face modes'); 

writeIntMatrix(fid, stackedNz', 'Sparse differentiation matrix ids'); 

writeFloatMatrix(fid, sparseSrr', 'Sparse differentiation Srr values');
writeFloatMatrix(fid, sparseSrs', 'Sparse differentiation Srs values');
writeFloatMatrix(fid, sparseSss', 'Sparse differentiation Sss values');

fclose(fid)

end
