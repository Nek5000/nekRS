function cub = CubatureVolumeMesh2D(CubatureOrder)

% function cub = CubatureVolumeMesh2D(CubatureOrder)
% purpose: build cubature nodes, weights and geometric factors for all elements

Globals2D;

% set up cubature nodes
[cub.R,cub.S,cub.W, cub.Ncub] = Cubature2D(CubatureOrder); 

% evaluate generalized Vandermonde of Lagrange interpolant functions at cubature nodes
cub.V  = InterpMatrix2D(cub.R, cub.S); 

% evaluate local derivatives of Lagrange interpolants at cubature nodes
[cub.Dr,cub.Ds] = Dmatrices2D(N,cub.R,cub.S,V);

% evaluate the geometric factors at the cubature nodes
[cub.rx,cub.sx,cub.ry,cub.sy,cub.J] = ...
    GeometricFactors2D(x,y, cub.Dr,cub.Ds);

% custom mass matrix per element
cub.mmCHOL = zeros(Np, Np, K);
cub.mm     = zeros(Np, Np, K);
for k=1:K
  cub.mm(:,:,k)     = cub.V'*diag(cub.J(:,k).*cub.W(:))*cub.V;
  cub.mmCHOL(:,:,k) = chol(cub.mm(:,:,k));
end

% incorporate weights and Jacobian
cub.w = cub.W; cub.W = cub.W*ones(1,K); cub.W = cub.W.*cub.J; 

% compute coordinates of cubature nodes
cub.x = cub.V*x;
cub.y = cub.V*y;

return
