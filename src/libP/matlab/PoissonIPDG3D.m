function [OP,MM] = PoissonIPDG3D()

% function [OP,MM] = PoissonIPDG3D()
% Purpose: Set up the discrete Poisson matrix directly
%          using IP. The operator is set up in the weak form

Globals3D;

% build local face matrices
massEdge = zeros(Np,Np,Nfaces);

% face mass matrix 1
Fm = Fmask(:,1); faceR = r(Fm); faceS = s(Fm); 
V2D = Vandermonde2D(N, faceR, faceS);
massEdge(Fm,Fm,1) = inv(V2D*V2D');

% face mass matrix 2
Fm = Fmask(:,2); faceR = r(Fm); faceT = t(Fm);
V2D = Vandermonde2D(N, faceR, faceT);
massEdge(Fm,Fm,2) = inv(V2D*V2D');

% face mass matrix 3
Fm = Fmask(:,3); faceS = s(Fm); faceT = t(Fm);
V2D = Vandermonde2D(N, faceS, faceT); 
massEdge(Fm,Fm,3) = inv(V2D*V2D');

% face mass matrix 4
Fm = Fmask(:,4); faceS = s(Fm); faceT = t(Fm);
V2D = Vandermonde2D(N, faceS, faceT); 
massEdge(Fm,Fm,4) = inv(V2D*V2D');

% build local volume mass matrix
MassMatrix = invV'*invV;

% build DG derivative matrices
MM  = zeros(K*Np*Np, 3);  OP = zeros(K*Np*Np*(1+Nfaces), 3);  

% global node numbering
entries = (1:Np*Np)'; entriesMM = (1:Np*Np)'; 
for k1=1:K 
  if(~mod(k1,1000)) k1, end;
  rows1 = ((k1-1)*Np+1:k1*Np)'*ones(1,Np); cols1 = rows1';

  % Build local operators  
  Dx = rx(1,k1)*Dr + sx(1,k1)*Ds + tx(1,k1)*Dt;   
  Dy = ry(1,k1)*Dr + sy(1,k1)*Ds + ty(1,k1)*Dt;
  Dz = rz(1,k1)*Dr + sz(1,k1)*Ds + tz(1,k1)*Dt;

  OP11 = J(1,k1)*(Dx'*MassMatrix*Dx + Dy'*MassMatrix*Dy + Dz'*MassMatrix*Dz);

  % Build element-to-element parts of operator
  for f1=1:Nfaces
    k2 = EToE(k1,f1); f2 = EToF(k1,f1); 

    rows2 = ((k2-1)*Np+1:k2*Np)'*ones(1,Np); cols2 = rows2';
    
    fidM  = (k1-1)*Nfp*Nfaces + (f1-1)*Nfp + (1:Nfp);
    vidM = vmapM(fidM); Fm1 = mod(vidM-1,Np)+1;
    vidP = vmapP(fidM); Fm2 = mod(vidP-1,Np)+1;
    
    id = 1+(f1-1)*Nfp + (k1-1)*Nfp*Nfaces;
    lnx = nx(id);  lny = ny(id);  lnz = nz(id); lsJ = sJ(id); 
    hinv = max(Fscale(id), Fscale(1+(f2-1)*Nfp, k2));    

    Dx2 = rx(1,k2)*Dr + sx(1,k2)*Ds + tx(1,k2)*Dt;   
    Dy2 = ry(1,k2)*Dr + sy(1,k2)*Ds + ty(1,k2)*Dt;
    Dz2 = rz(1,k2)*Dr + sz(1,k2)*Ds + tz(1,k2)*Dt;
    
    Dn1 = lnx*Dx  + lny*Dy  + lnz*Dz;
    Dn2 = lnx*Dx2 + lny*Dy2 + lnz*Dz2;

    mmE = lsJ*massEdge(:,:,f1);

    gtau = 2*(N+1)*(N+1)*hinv; % set penalty scaling
    
    if(EToE(k1,f1)==k1)
      OP11 = OP11 + ( gtau*mmE - mmE*Dn1 - Dn1'*mmE ); % ok
    else 
      % interior face variational terms
      OP11        = OP11 + 0.5*( gtau*mmE - mmE*Dn1 - Dn1'*mmE );
      
      OP12 = zeros(Np);
      OP12(:,Fm2) =             - 0.5*( gtau*mmE(:,Fm1) );
      OP12(Fm1,:) = OP12(Fm1,:) - 0.5*(      mmE(Fm1,Fm1)*Dn2(Fm2,:) );
      OP12(:,Fm2) = OP12(:,Fm2) - 0.5*(-Dn1'*mmE(:, Fm1) );
      OP(entries(:), :) = [rows1(:), cols2(:), OP12(:)];
      entries = entries + Np*Np;
    end 
  end      
  OP(entries(:), :)   = [rows1(:), cols1(:), OP11(:)];
  MM(entriesMM(:), :) = [rows1(:), cols1(:), J(1,k1)*MassMatrix(:)];
  entries = entries + Np*Np; entriesMM = entriesMM + Np*Np;
end  

OP   =   OP(1:max(entries)  -Np*Np,:);  OP   = myspconvert(OP, Np*K, Np*K, 1e-15);
MM   =   MM(1:max(entriesMM)-Np*Np,:);  MM   = myspconvert(MM, Np*K, Np*K, 1e-15);
return
