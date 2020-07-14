function [vmapM, vmapP, vmapB, mapB] = BuildMaps3D

% function [vmapM, vmapP, vmapB, mapB] = BuildMaps3D
% Purpose: Connectivity and boundary tables for nodes given in the K # of elements,
% 	       each with N+1 degrees of freedom.

Globals3D;

% number volume nodes consecutively
nodeids = reshape(1:K*Np, Np, K);
vmapM   = zeros(Nfp, Nfaces, K); 
vmapP   = zeros(Nfp, Nfaces, K); 

for k1=1:K
  for f1=1:Nfaces
    % find index of face nodes with respect to volume node ordering
    vmapM(:,f1,k1) = nodeids(Fmask(:,f1), k1);
  end
end

tmp = ones(1,Nfp);
for k1=1:K
  for f1=1:Nfaces
    % find neighbor
    k2 = EToE(k1,f1); f2 = EToF(k1,f1);
    
    % find find volume node numbers of left and right nodes 
    vidM = vmapM(:,f1,k1); vidP = vmapM(:,f2,k2);
    
    xM = x(vidM)*tmp; yM = y(vidM)*tmp; zM = z(vidM)*tmp;
    xP = x(vidP)*tmp; yP = y(vidP)*tmp; zP = z(vidP)*tmp;
    
    % Compute distance matrix
    D = (xM -xP').^2 + (yM-yP').^2 + (zM-zP').^2;

    [idM, idP] = find(abs(D)<NODETOL);
    vmapP(idM, f1, k1) = vmapM(idP, f2, k2);
  end
end

vmapP = vmapP(:); vmapM = vmapM(:);

% Create list of boundary nodes
mapB = find(vmapP==vmapM); vmapB = vmapM(mapB);
return
