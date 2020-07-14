function [Np, EToV, VX, V, Dr, x, drdx, invM, vmapM, vmapP, mapB, nx, LIFT] = ...
  startUp1D(N, Nelements, xmin, xmax)

Np = N+1;

%% set up nodes 
r = JacobiGL(0,0,N);    %% GLL nodes

%% operator matrices
V = Vandermonde1D(N, r); %% Legendre Vandermonde
Dr = Dmatrix1D(N, r, V); %% Collocation differentiation matrix
M = inv(V')/V;
M = diag(sum(M));
invM = inv(M);
%invM = V*transpose(V);   %% Inverse mass matrix

l0 = [1;zeros(N,1)];     %% Left Lagrange node
lN = [zeros(N,1);1];     %% Right Lagrange node
LIFT = invM*[l0,lN];     %% LIFT matrix

%% mesh
Nnodes = Nelements+1;
VX = linspace(xmin, xmax, Nnodes);
EToV = [(1:Nelements)', (2:Nelements+1)'];

%% geometric factors 
x = 0.5*(1-r)*VX(EToV(:,1)) + 0.5*(1+r)*VX(EToV(:,2));
dxdr = Dr*x;
drdx = 1./dxdr;
nx = [-1;1]*ones(1,Nelements);

%% arrays of indices for trace nodes 
ids = reshape(1:Np*Nelements, Np, Nelements);
vmapM = [ids(1,:);ids(Np,:)];
vmapP = vmapM;
for e=1:Nelements
  if(e>1)  %% face 1
    vmapP(1,e) = ids(Np,e-1); 
  end;
  if(e<Nelements) %% face 2
    vmapP(2,e) = ids(1,e+1); 
  end;
end

%% boundary trace nodes
mapB = [1,2*Nelements];
