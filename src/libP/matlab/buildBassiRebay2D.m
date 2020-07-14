Globals2D

%% degree
N = 2;

%% Read in Mesh
[Nv, VX, VY, K, EToV] = MeshReaderGmsh2D('../meshes/tri.msh');

%% Initialize solver and construct grid and metric
StartUp2D;

%% default to Dirichlet
BCType = Dirichlet*(EToE == (1:K)'*ones(1,Nfaces));

%% set up boundary conditions
BuildBCMaps2D;

%% build reference mass matrix
MassMatrix = inv(V*V');

%% build matrices
dxD = zeros(Np*K);
dyD = zeros(Np*K);
dxN = zeros(Np*K);
dyN = zeros(Np*K);
MM  = zeros(Np*K);
S   = zeros(Np*K);

for k=1:K

  dx = rx(1,k)*Dr + sx(1,k)*Ds;
  dy = ry(1,k)*Dr + sy(1,k)*Ds;

  ids = (1:Np) + (k-1)*Np;

  dxD(ids,ids) = dxD(ids,ids) + dx;
  dyD(ids,ids) = dyD(ids,ids) + dy;

  dxN(ids,ids) = dxN(ids,ids) + dx;
  dyN(ids,ids) = dyN(ids,ids) + dy;

  MM(ids,ids) = J(1,k)*MassMatrix;

  for f=1:Nfaces
    nids = (1+(f-1)*Nfp):f*Nfp;
    fids = (k-1)*Nfp*Nfaces + nids;
    idsM = vmapM(fids);
    idsP = vmapP(fids);

    bc = BCType(k,f);

    LIFTx = LIFT(:,nids)*diag(nx(fids).*sJ(fids)/J(1,k));
    LIFTy = LIFT(:,nids)*diag(ny(fids).*sJ(fids)/J(1,k));
    STAB  = LIFT(:,nids)*diag(sJ(fids)/J(1,k));
    switch bc
      case Dirichlet
  dxD(ids,idsM) = dxD(ids,idsM) - LIFTx;
  dyD(ids,idsM) = dyD(ids,idsM) - LIFTy;

  S(ids,idsM) = S(ids,idsM)-STAB;
      case Neuman
  dxN(ids,idsM) = dxN(ids,idsM) - LIFTx;
  dyN(ids,idsM) = dyN(ids,idsM) - LIFTy;

      otherwise
  dxD(ids,idsM) = dxD(ids,idsM) - 0.5*LIFTx;
  dyD(ids,idsM) = dyD(ids,idsM) - 0.5*LIFTy;

  dxD(ids,idsP) = dxD(ids,idsP) + 0.5*LIFTx;
  dyD(ids,idsP) = dyD(ids,idsP) + 0.5*LIFTy;

  dxN(ids,idsM) = dxN(ids,idsM) - 0.5*LIFTx;
  dyN(ids,idsM) = dyN(ids,idsM) - 0.5*LIFTy;

  dxN(ids,idsP) = dxN(ids,idsP) + 0.5*LIFTx;
  dyN(ids,idsP) = dyN(ids,idsP) + 0.5*LIFTy;

  S(ids,idsM) = S(ids,idsM)-0.5*STAB;
  S(ids,idsP) = S(ids,idsP)+0.5*STAB;

    end
  end
end

%% this is symmetric, positive definite
tau = 1;
A = -MM*(dxN*dxD + dyN*dyD + tau*S);