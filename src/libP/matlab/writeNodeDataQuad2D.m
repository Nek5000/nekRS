
function writeNodeDataQuad2D(N)

r1d = JacobiGL(0,0,N);
V1d = Vandermonde1D(N, r1d);
D1d = Dmatrix1D(N, r1d, V1d);
M1d = inv(V1d')/V1d;
w1d = sum(M1d);
cnt = 1;
for i=1:N+1
  for j=1:N+1
    r(cnt) = r1d(j); %% r runs fastest
    s(cnt) = r1d(i);
    cnt = cnt+1;
  end
end
Np = (N+1)*(N+1);
r = reshape(r, Np,1);
s = reshape(s, Np,1);

Nfp = N+1;
Nfaces = 4;

% find all the nodes that lie on each edge
NODETOL = 1e-8;
faceNodes1   = find( abs(s+1) < NODETOL)'; 
faceNodes2   = find( abs(r-1) < NODETOL)';
faceNodes3   = find( abs(s-1) < NODETOL)';
faceNodes4   = find( abs(r+1) < NODETOL)';
faceNodes  = [faceNodes1;faceNodes2;faceNodes3;faceNodes4]';

V = VandermondeQuad2D(N, r, s);
[Dr,Ds] = DmatricesQuad2D(N, r, s, V);
Dr = full(Dr);
Ds = full(Ds);
LIFT = LiftQuad2D(N, faceNodes, r, s);
LIFT = full(LIFT);

fname = sprintf('quadrilateralN%02d.dat', N);
fid = fopen(fname, 'w');

writeFloatMatrix(fid, r, 'Nodal r-coordinates');
writeFloatMatrix(fid, s, 'Nodal s-coordinates');
writeFloatMatrix(fid, Dr, 'Nodal Dr differentiation matrix');
writeFloatMatrix(fid, Ds, 'Nodal Ds differentiation matrix');

writeIntMatrix(fid, faceNodes'-1, 'Nodal Face nodes');
writeFloatMatrix(fid, LIFT, 'Nodal Lift Matrix');

writeFloatMatrix(fid, r1d, 'Nodal 1D GLL Nodes');
writeFloatMatrix(fid, w1d', 'Nodal 1D GLL Weights');
writeFloatMatrix(fid, D1d, 'Nodal 1D differentiation matrix');


%% compute equispaced nodes on equilateral triangle
[plotR,plotS] = meshgrid(linspace(-1,1,N+4));
plotR = plotR(:); plotS = plotS(:);

%% count plot nodes
plotNp = length(plotR);

%% triangulate equilateral element nodes
plotEToV = delaunay(plotR,plotS)-1; 

%% count triangles in plot node triangulation
plotNelements = size(plotEToV,1); 

%check triangulation
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

%% create interpolation matrix from warp & blend to plot nodes
plotInterp = VandermondeQuad2D(N, plotR,plotS)/V; 

writeFloatMatrix(fid, plotR, 'Plotting r-coordinates');
writeFloatMatrix(fid, plotS, 'Plotting s-coordinates');
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

%% 1D quadrature
[z,w] = JacobiGQ(0,0,ceil(3*N/2));
cV1d = Vandermonde1D(N, z);
cVr = GradVandermonde1D(N, z);
cInterp = cV1d/V1d;
cDiff = cVr/V1d;
cubProject = (cV1d/V1d)'*diag(w);
cubDT = (cVr/V1d)'*diag(w);
Nqc = length(z);

writeFloatMatrix(fid, z, 'Quadrature r-coordinates');
writeFloatMatrix(fid, w, 'Quadrature weights');

writeFloatMatrix(fid, cInterp, 'Quadrature Interpolation Matrix');
writeFloatMatrix(fid, cDiff,   'Quadrature Differentiation Interpolation Matrix');
writeFloatMatrix(fid, cubDT, 'Quadrature Weak D Differentiation Matrix');
writeFloatMatrix(fid, cubProject, 'Quadrature Projection Matrix');


%{
%% volume cubature
[z,w] = JacobiGQ(0,0,ceil(3*N/2));
[cubr,cubs] = meshgrid(z);
cubw = w*transpose(w);
cubr = cubr(:);
cubs = cubs(:);
cubw = cubw(:);

cInterp = VandermondeQuad2D(N, cubr, cubs)/V;
Ncub = length(cubr);

cV = VandermondeQuad2D(N, cubr, cubs);
cV'*diag(cubw)*cV;

[cVr,cVs] = GradVandermondeQuad2D(N, cubr, cubs);
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


%% surface cubature
[z,w] = JacobiGQ(0,0,ceil(3*N/2));
%z = JacobiGL(0,0,N);
%zV = Vandermonde1D(N,z);
%w = sum(inv(zV*transpose(zV)));

Nfi = length(z);

ir = [z,ones(Nfi,1),-z,-ones(Nfi,1)];
is = [-ones(Nfi,1), z, ones(Nfi,1), -z];
iw = [w,w,w,w];

sV = VandermondeQuad2D(N, ir(:), is(:));
sInterp = sV/V;

interp = [sInterp(1:Nfi,faceNodes(:,1));
sInterp(Nfi+1:2*Nfi,faceNodes(:,2));
sInterp(2*Nfi+1:3*Nfi,faceNodes(:,3));
sInterp(3*Nfi+1:4*Nfi,faceNodes(:,4))];

iLIFT = V*V'*sInterp'*diag(iw(:));

writeFloatMatrix(fid, interp, 'Cubature Surface Interpolation Matrix');
writeFloatMatrix(fid, iLIFT, 'Cubature Surface Lift Matrix');
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
FEMEToV = zeros(N*N,4);
cnt =1;
for j=0:N-1
  for i=0:N-1
    FEMEToV(cnt,:) = [i+(j)*(N+1), i+1+(j)*(N+1), i+1+(j+1)*(N+1), i+(j+1)*(N+1)];
    cnt = cnt+1;
  end
end

writeIntMatrix(fid, FEMEToV, 'SEMFEM reference mesh');

fclose(fid);
