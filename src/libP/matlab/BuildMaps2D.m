function [mapM, mapP, vmapM, vmapP, vmapB, mapB] = BuildMaps2D()

% function [mapM, mapP, vmapM, vmapP, vmapB, mapB] = BuildMaps2D
% Purpose: Connectivity and boundary tables for nodes given
% 	   in the K # of elements, each with Np degrees of freedom.

Globals2D;

% number volume nodes consecutively
nodeids = reshape(1:K*Np, Np, K);
vmapM   = zeros(Nfp, Nfaces, K); 
vmapP   = zeros(Nfp, Nfaces, K); 
mapM    = (1:K*Nfp*Nfaces)';
mapP    = reshape(mapM, Nfp, Nfaces, K);
 
% find index of face nodes with respect to volume node ordering
for k1=1:K
  for f1=1:Nfaces
    vmapM(:,f1,k1) = nodeids(Fmask(:,f1), k1);
  end
end

one = ones(1, Nfp);
for k1=1:K
  for f1=1:Nfaces
    % find neighbor
    k2 = EToE(k1,f1); f2 = EToF(k1,f1);
    
    % reference length of edge
    v1 = EToV(k1,f1); v2 = EToV(k1, 1+mod(f1,Nfaces));
    refd = sqrt( (VX(v1)-VX(v2))^2 + (VY(v1)-VY(v2))^2 );

    % find find volume node numbers of left and right nodes 
    vidM = vmapM(:,f1,k1);
    vidP = vmapM(:,f2,k2);
    
    x1 = x(vidM); y1 = y(vidM); x2 = x(vidP); y2 = y(vidP);
    x1 = x1*one;  y1 = y1*one;  x2 = x2*one;  y2 = y2*one;

    % Compute distance matrix
    D = (x1 -x2').^2 + (y1-y2').^2;
    [idM, idP] = find(sqrt(abs(D))<NODETOL*refd);
    vmapP(idM,f1,k1) = vidP(idP); 
    mapP(idM,f1,k1) = idP + (f2-1)*Nfp+(k2-1)*Nfaces*Nfp;
  end
end

% reshape vmapM and vmapP to be vectors
vmapP = vmapP(:); vmapM = vmapM(:); mapP = mapP(:); 
% Create list of boundary nodes
mapB = find(vmapP==vmapM); vmapB = vmapM(mapB);
return
