% Driver script for solving the 2D vacuum Maxwell's equations on TM form
Globals2D;

% Polynomial order used for approximation 
N = 1;

% Read in Mesh
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell025.neu');
%[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('maxmesh.neu');

[Nv, VX, VY, K, EToV] = MeshReaderGmsh2D('cavityH01.msh');

[r,s] = Nodes2D(N);
[r,s] = xytors(r,s);

Np = length(r);
Nfp = N+1;
Nfaces = 3;

% find all the nodes that lie on each edge
NODETOL = 1e-8;
faceNodes1   = find( abs(s+1) < NODETOL)'; 
faceNodes2   = find( abs(r+s) < NODETOL)';
faceNodes3   = find( abs(r+1) < NODETOL)';
FaceNodes  = [faceNodes1;faceNodes2;faceNodes3]';
Fmask = FaceNodes;

% build coordinates of all the nodes
va = EToV(:,1)'; vb = EToV(:,2)'; vc = EToV(:,3)';
x = 0.5*(-(r+s)*VX(va)+(1+r)*VX(vb)+(1+s)*VX(vc));
y = 0.5*(-(r+s)*VY(va)+(1+r)*VY(vb)+(1+s)*VY(vc));

V = Vandermonde2D(N, r, s);
[Dr,Ds] = Dmatrices2D(N, r, s, V);
LIFT = Lift2D(N, FaceNodes, r, s);

% calculate geometric factors
[rx,sx,ry,sy,J] = GeometricFactors2D(x,y,Dr,Ds);

% calculate geometric factors
[nx, ny, sJ] = Normals2D();
Fscale = sJ./(J(Fmask,:));


%% volume cubature
[cubr,cubs,cubw] = Cubature2D(3*N);
cInterp = Vandermonde2D(N, cubr, cubs)/V;
Ncub = length(cubr);

cV = Vandermonde2D(N, cubr, cubs);
cV'*diag(cubw)*cV

[cVr,cVs] = GradVandermonde2D(N, cubr, cubs);
cubDrT = V*transpose(cVr)*diag(cubw);
cubDsT = V*transpose(cVs)*diag(cubw);

%% surface cubature
[z,w] = JacobiGQ(0,0,ceil(3*N/2));
%z = JacobiGL(0,0,N);
%zV = Vandermonde1D(N,z);
%w = sum(inv(zV*transpose(zV)));

Nfi = length(z);

ir = [z,-z,-ones(Nfi,1)];
is = [-ones(Nfi,1), z, -z];
iw = [w,w,w];

sV = Vandermonde2D(N, ir(:), is(:));
sInterp = sV/V;

interp = [sInterp(1:Nfi,FaceNodes(:,1));sInterp(Nfi+1:2*Nfi,FaceNodes(:,2));sInterp(2*Nfi+1:3*Nfi,FaceNodes(:,3))];

iLIFT = V*V'*sInterp'*diag(iw(:));

nx = [ones(Nfi,1)*nx(1,:);ones(Nfi,1)*nx(Nfp+1,:);ones(Nfi,1)*nx(2*Nfp+1,:)];
ny = [ones(Nfi,1)*ny(1,:);ones(Nfi,1)*ny(Nfp+1,:);ones(Nfi,1)*ny(2*Nfp+1,:)];
sJ = [ones(Nfi,1)*sJ(1,:);ones(Nfi,1)*sJ(Nfp+1,:);ones(Nfi,1)*sJ(2*Nfp+1,:)];
rx = [ones(Ncub,1)*rx(1,:)];
sx = [ones(Ncub,1)*sx(1,:)];
ry = [ones(Ncub,1)*ry(1,:)];
sy = [ones(Ncub,1)*sy(1,:)];
J  = [ones(Ncub,1)*J(1,:)];
iJ = [ones(Nfi*Nfaces,1)*J(1,:)];
vcon = cInterp*ones(Np,K);
scon = sInterp*ones(Np,K); % ones(Nfi*Nfaces,K);
cubDrT*(rx.*vcon)+cubDsT*(sx.*vcon) - iLIFT*((sJ./iJ).*nx.*scon)
cubDrT*(ry.*vcon)+cubDsT*(sy.*vcon) - iLIFT*((sJ./iJ).*ny.*scon)

