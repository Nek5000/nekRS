
function writeNodeData3D(inN)

Globals3D;
N = inN;
[r,s,t] = Nodes3D(N);
[r,s,t] = xyztorst(r,s,t);

Np = length(r);
Nfp = (N+1)*(N+2)/2;
Nfaces = 4;

%% Nodal data
% find all the nodes that lie on each edge
NODETOL = 1e-8;
faceNodes1   = find( abs(t+1) < NODETOL)';
faceNodes2   = find( abs(s+1) < NODETOL)';
faceNodes3   = find( abs(r+s+t+1) < NODETOL)';
faceNodes4   = find( abs(r+1) < NODETOL)';
faceNodes  = [faceNodes1;faceNodes2;faceNodes3;faceNodes4]';

V = Vandermonde3D(N, r, s, t);
[Dr,Ds,Dt] = Dmatrices3D(N, r, s, t, V);
MM = inv(transpose(V))/V;
LIFT = Lift3D(N, faceNodes, r, s, t);

fname = sprintf('tetN%02d.dat', N);
fid = fopen(fname, 'w');

writeFloatMatrix(fid, r, 'Nodal r-coordinates');
writeFloatMatrix(fid, s, 'Nodal s-coordinates');
writeFloatMatrix(fid, t, 'Nodal t-coordinates');
writeFloatMatrix(fid, Dr, 'Nodal Dr differentiation matrix');
writeFloatMatrix(fid, Ds, 'Nodal Ds differentiation matrix');
writeFloatMatrix(fid, Dt, 'Nodal Dt differentiation matrix');
writeFloatMatrix(fid, MM, 'Nodal Mass Matrix');

writeIntMatrix(fid, faceNodes'-1, 'Nodal Face nodes');
writeFloatMatrix(fid, LIFT, 'Nodal Lift Matrix');


%% Plotting data
% compute equispaced nodes on equilateral triangle
[plotR,plotS,plotT] = EquiNodes3D(N+4);
plotNp = length(plotR);
plotEToV = delaunayFixVolume(plotR,plotS,plotT)-1;
plotNelements = size(plotEToV,1);
plotInterp = Vandermonde3D(N, plotR,plotS,plotT)/V;

writeFloatMatrix(fid, plotR, 'Plotting r-coordinates');
writeFloatMatrix(fid, plotS, 'Plotting s-coordinates');
writeFloatMatrix(fid, plotT, 'Plotting t-coordinates');
writeFloatMatrix(fid, plotInterp, 'Plotting Interpolation Matrix');
writeIntMatrix(fid, plotEToV, 'Plotting triangulation');

%% Cubature data 
if N < 7
    [cubr cubs cubt cubw] = tet_cubature(2*N+1);
    Vq = Vandermonde3D(N,cubr,cubs,cubt)/V;
    Pq = V*V'*Vq'*diag(cubw);
    
    cV = Vandermonde3D(N, cubr, cubs, cubt);
    [cVr,cVs,cVt] = GradVandermonde3D(N, cubr, cubs, cubt);
    cubDrT = V*transpose(cVr)*diag(cubw);
    cubDsT = V*transpose(cVs)*diag(cubw);
    cubDtT = V*transpose(cVt)*diag(cubw);
    
    writeFloatMatrix(fid, cubr, 'Cubature r-coordinates');
    writeFloatMatrix(fid, cubs, 'Cubature s-coordinates');
    writeFloatMatrix(fid, cubt, 'Cubature t-coordinates');
    writeFloatMatrix(fid, cubw, 'Cubature weights');

    writeFloatMatrix(fid, Vq, 'Cubature Interpolation Matrix');
    writeFloatMatrix(fid, cubDrT, 'Cubature Weak Dr Differentiation Matrix');
    writeFloatMatrix(fid, cubDsT, 'Cubature Weak Ds Differentiation Matrix');
    writeFloatMatrix(fid, cubDtT, 'Cubature Weak Dt Differentiation Matrix');
    writeFloatMatrix(fid, Pq, 'Cubature Projection Matrix');  
    
    % Surface Cubature
    [cubx,cuby,cubw] = Cubature2D(3*N); 
    Nfi = length(cubx);
    ir = [cubx,         cubx,        cubx,                  -ones(Nfi,1)];
    is = [cuby,        -ones(Nfi,1), cuby,                   cubx];
    it = [-ones(Nfi,1), cuby,       -(ones(Nfi,1)+cubx+cuby),cuby];
    iw = [cubw,cubw,cubw,cubw];
    
    sV = Vandermonde3D(N, ir(:), is(:),it(:));
    sInterp = sV/V;
    
    interp = [sInterp(1:Nfi,faceNodes(:,1));...
        sInterp(Nfi+1:2*Nfi,faceNodes(:,2));...
        sInterp(2*Nfi+1:3*Nfi,faceNodes(:,3));...
        sInterp(3*Nfi+1:4*Nfi,faceNodes(:,4))];
    
    % integration node lift matrix
    iLIFT = V*V'*sInterp'*diag(iw(:));
    
    writeFloatMatrix(fid, interp, 'Cubature Surface Interpolation Matrix');
    writeFloatMatrix(fid, iLIFT, 'Cubature Surface Lift Matrix');
    
end

%% Berstein-Bezier operators

addpath('./bern')
tol=1e-6;

[VB,Vr,Vs,Vt,V0,V1,V2,V3] = bern_basis_tet(N,r,s,t);

invVB = inv(VB);

[D0ids,D0vals,D1ids,D1vals,D2ids,D2vals,D3ids,D3vals] = bern_basis_diff3D(N,V,VB,V0,V1,V2,V3);
[D0Tids,D0Tvals,D1Tids,D1Tvals,D2Tids,D2Tvals,D3Tids,D3Tvals] = bern_basis_transposediff3D(N,V,VB,V0,V1,V2,V3);
[L0ids,L0vals,ELids,ELvals] = bern_basis_lift3D(N,V,VB,r,s,t);

%write out the BB operators
writeFloatMatrix(fid, VB, 'Bernstein-Bezier Vandermonde Matrix');
writeFloatMatrix(fid, invVB, 'Bernstein-Bezier Inverse Vandermonde Matrix');
writeIntMatrix(fid, D0ids, 'Bernstein-Bezier sparse D0 differentiation ids');
writeIntMatrix(fid, D1ids, 'Bernstein-Bezier sparse D1 differentiation ids');
writeIntMatrix(fid, D2ids, 'Bernstein-Bezier sparse D2 differentiation ids');
writeIntMatrix(fid, D3ids, 'Bernstein-Bezier sparse D3 differentiation ids');
writeFloatMatrix(fid, D0vals, 'Bernstein-Bezier sparse D differentiation values');

writeIntMatrix(fid, D0Tids, 'Bernstein-Bezier sparse D0T transpose differentiation ids');
writeIntMatrix(fid, D1Tids, 'Bernstein-Bezier sparse D1T transpose differentiation ids');
writeIntMatrix(fid, D2Tids, 'Bernstein-Bezier sparse D2T transpose differentiation ids');
writeIntMatrix(fid, D3Tids, 'Bernstein-Bezier sparse D3T transpose differentiation ids');
writeFloatMatrix(fid, D0Tvals, 'Bernstein-Bezier sparse DT transpose differentiation values');

writeIntMatrix(fid, L0ids, 'Bernstein-Bezier L0 Matrix ids');
writeFloatMatrix(fid, L0vals, 'Bernstein-Bezier L0 Matrix values');
writeIntMatrix(fid, ELids, 'Bernstein-Bezier EL lift ids');
writeFloatMatrix(fid, ELvals, 'Bernstein-Bezier EL lift values');

%degree raise/lower operators along traces
[r2Dp1 s2Dp1] = Nodes2D(N+1); [r2Dp1 s2Dp1] = xytors(r2Dp1,s2Dp1);
BBRaise = bern_basis_tri(N+1,r2Dp1,s2Dp1)\bern_basis_tri(N,r2Dp1,s2Dp1);
BBRaise(abs(BBRaise)<tol) = 0;

Nfpp1 = (N+2)*(N+3)/2;
BBRaiseIds  = zeros(Nfpp1,3);
BBRaiseVals = zeros(Nfpp1,3);

for i = 1:Nfpp1
    tmp = find(BBRaise(i,:));
    BBRaiseVals(i,1:length(tmp)) = BBRaise(i,tmp);
    tmp = tmp-1; % zero indexing
    if length(tmp) < 3
        tmp = [tmp zeros(1,3-length(tmp))];
    end
    BBRaiseIds(i,:) = tmp;
end

[r2D,s2D] = Nodes2D(N); [r2D,s2D] = xytors(r2D,s2D);
VB2D = bern_basis_tri(N,r2D,s2D);
V2D = Vandermonde2D(N, r2D,s2D);

Nm1 = N-1;
if (N>1)
    [r2Dm1 s2Dm1] = Nodes2D(Nm1); [r2Dm1 s2Dm1] = xytors(r2Dm1,s2Dm1);
    VB2Dm1 = bern_basis_tri(Nm1,r2Dm1,s2Dm1);
    V2Dm1 = Vandermonde2D(Nm1, r2Dm1,s2Dm1);
else
    r2Dm1 =0;
    s2Dm1 =0;
    VB2Dm1 = 1;
    V2Dm1 = 1;
end

Nfpm1 = (Nm1+1)*(Nm1+2)/2;

BBLower = V2D\VB2D;
BB = [];
sk =1;
for n=0:N
    for m=0:N-n
        if n+m<N
            BB = [BB; BBLower(sk,:)];
        end
        sk = sk+1;
    end
end
BBLower = VB2Dm1\V2Dm1*BB;

writeIntMatrix(fid, BBRaiseIds, 'Bernstein-Bezier sparse 2D degree raise ids');
writeFloatMatrix(fid, BBRaiseVals, 'Bernstein-Bezier sparse 2D degree raise values');
writeFloatMatrix(fid, BBLower, 'Bernstein-Bezier sparse 2D degree lower matrix');


%% elliptic patch problem
K = 5;

VX = [-1, 1,      0,          0,           0,         5/3,        -5/3,         0];
VY = [ 0, 0,sqrt(3),  1/sqrt(3),-7*sqrt(3)/9, 8*sqrt(3)/9, 8*sqrt(3)/9, 1/sqrt(3)];
VZ = [ 0, 0,      0,2*sqrt(6)/3, 4*sqrt(6)/9, 4*sqrt(6)/9, 4*sqrt(6)/9,-2*sqrt(6)/3];

EToV = [1,2,3,4;
    1,2,4,5;
    2,3,4,6;
    3,1,4,7;
    1,3,2,8];

BCType = [0,0,0,0;
    0,0,0,0;
    0,0,0,0;
    0,0,0,0;
    0,0,0,0];

StartUp3D;

% build weak Poisson operator matrices
[A, M] = PoissonIPDG3D();

%% hack since we know W&B face 1 nodes are first
vmapP = reshape(vmapP, Nfp*Nfaces, K);
idsP = vmapP(:,1);
subind = [(1:Np)';idsP];

subA = full(A(subind,subind));
subM = full(M(subind,subind));

%spy(abs(subA)>1e-10);
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


%% permutations
mVXYZ = [-1,+1,-1,-1;
	 -1,-1,+1,-1;
	 -1,-1,-1,+1];
Nverts = 4;
cnt = 0;
for v1=1:Nverts
  for v2=1:Nverts
    for v3=1:Nverts
      for v4=1:Nverts

	if(v1~=v2 & v1~=v3 & v1~=v4 & v2~=v3 & v2 ~=v4 & v3~=v4)
	  vX1 = mVXYZ(:,v1)';
	  vX2 = mVXYZ(:,v2)';
	  vX3 = mVXYZ(:,v3)';
	  vX4 = mVXYZ(:,v4)';
				  
	  permRST = -0.5*(1+r+s+t)*vX1+0.5*(1+r)*vX2+0.5*(1+s)*vX3+0.5*(1+t)*vX4;

	  %[v1,v2,v3,v4]
	  
	  permr = permRST(:,1);
	  perms = permRST(:,2);
	  permt = permRST(:,3);

	  for n=1:Np
	    for m=1:Np
	      dist(n,m) = (r(n)-permr(m))^2 + (s(n)-perms(m))^2 + (t(n)-permt(m))^2;
	    end
	  end
	  [foo,ids] = min(dist);
	  cnt = cnt+1;
	  %ids
	  pmap(cnt, :) = [v1,v2,v3,v4,ids];

	end
      end
    end
  end
end


%% degree raising interpolation
[rP1,sP1,tP1] = Nodes3D(N+1);
[rP1,sP1,tP1] = xyztorst(rP1,sP1,tP1);

VP1 = Vandermonde3D(N, rP1, sP1, tP1);
IP1 = VP1/V;
NpP1 = length(rP1);

%% degree lowering interpolation
if(N>1)
  [rM1,sM1,tM1] = Nodes3D(N-1);
[rM1,sM1,tM1] = xyztorst(rM1,sM1,tM1);
else
%% hard code degree 0
rM1 = -1/2; % -1-1-1+1/4
sM1 = -1/2;
tM1 = -1/2;
end

VM1 = Vandermonde3D(N, rM1, sM1, tM1);
IM1 = VM1/V;
NpM1 = length(rM1);


writeFloatMatrix(fid, IP1, 'Nodal degree raise matrix');
writeFloatMatrix(fid, IM1, 'Nodal degree lower matrix');

%% SEMFEM
[req,seq,teq] = EquiNodes3D(N);
FEMEToV = delaunayFixVolume(req,seq,teq)-1;;
rFEM = r;
sFEM = s;
tFEM = t;

NpFEM = length(rFEM);
NelFEM = size(FEMEToV,1);

IQN = Vandermonde3D(N, rFEM, sFEM, tFEM)/V;
invIQN = (transpose(IQN)*IQN)\(transpose(IQN));

writeFloatMatrix(fid, rFEM, 'SEMFEM r-coordinates');
writeFloatMatrix(fid, sFEM, 'SEMFEM s-coordinates');
writeFloatMatrix(fid, tFEM, 'SEMFEM t-coordinates');

writeIntMatrix(fid, FEMEToV, 'SEMFEM reference mesh');  
writeFloatMatrix(fid, invIQN', 'SEMFEM interpolation matrix');


% build interpolation matrix (coarse->fine)
EToVi = [1 5 7 8; 5 2 6 9; 7 6 3 10; 8 9 10 4; 8 5 7 9; 7 5 6 9; 8 9 7 10; 9 6 7 10];
VXi   = [-1  1 -1 -1  0  0 -1 -1  0 -1];
VYi   = [-1 -1  1 -1 -1  0  0 -1 -1  0];
VZi   = [-1 -1 -1  1 -1 -1 -1  0  0  0];

v1 = EToVi(:,1); v2 = EToVi(:,2); v3 = EToVi(:,3); v4 = EToVi(:,4);
ri = 0.5*(-(r+s+t+1)*VXi(v1) + (1+r)*VXi(v2) + (1+s)*VXi(v3) + (1+t)*VXi(v4) );
si = 0.5*(-(r+s+t+1)*VYi(v1) + (1+r)*VYi(v2) + (1+s)*VYi(v3) + (1+t)*VYi(v4) );
ti = 0.5*(-(r+s+t+1)*VZi(v1) + (1+r)*VZi(v2) + (1+s)*VZi(v3) + (1+t)*VZi(v4) );

contourInterp = Vandermonde3D(N, ri(:), si(:), ti(:))*invV;
ri = [-1;1;-1;-1]; si = [-1;-1;1;-1]; ti = [-1;-1;-1;1]; refNp = length(ri);
contourInterp1 = Vandermonde3D(N, ri(:), si(:), ti(:))*invV;

sk = 1;
F = zeros(Np);
for i=0:N % old ordering
  for j=0:N - i
    for k=0:N - i - j
      if(i+j+k<=1), F(sk,sk) = 1.; end;
      sk = sk+1;
    end
  end
end

contourFilter = V*F*invV;

writeIntMatrix(fid, EToVi, 'Contour plot EToV');
writeFloatMatrix(fid, VXi, 'Contour plot VX');
writeFloatMatrix(fid, VYi, 'Contour plot VY');
writeFloatMatrix(fid, VZi, 'Contour plot VZ');

writeFloatMatrix(fid, contourInterp, 'Contour plot Interpolation');
writeFloatMatrix(fid, contourInterp1, 'Contour plot Linear Interpolation');
writeFloatMatrix(fid, contourFilter, 'Contour plot Filter');

fclose(fid);
