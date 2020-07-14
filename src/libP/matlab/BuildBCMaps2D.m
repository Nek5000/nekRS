function BuildBCMaps2D()

% function BuildMaps2DBC
% Purpose: Build specialized nodal maps for various types of
%          boundary conditions, specified in BCType.

Globals2D;

% create label of face nodes with boundary types from BCType
bct    = BCType';
bnodes = ones(Nfp, 1)*bct(:)';
bnodes = bnodes(:);

% find location of boundary nodes in face and volume node lists
mapI = find(bnodes==In);           vmapI = vmapM(mapI);
mapO = find(bnodes==Out);          vmapO = vmapM(mapO);
mapW = find(bnodes==Wall);         vmapW = vmapM(mapW);
mapF = find(bnodes==Far);          vmapF = vmapM(mapF);
mapC = find(bnodes==Cyl);          vmapC = vmapM(mapC);
mapD = find(bnodes==Dirichlet);    vmapD = vmapM(mapD);
mapN = find(bnodes==Neuman);       vmapN = vmapM(mapN);
mapS = find(bnodes==Slip);         vmapS = vmapM(mapS);
return;