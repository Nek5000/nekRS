function [EToE,EToF]= tiConnect2D(EToV)

% function [EToE,EToF]= tiConnect2D(EToV)
% Purpose: triangle face connect algorithm due to Toby Isaac

Nfaces=3;
K = size(EToV,1);
Nnodes = max(max(EToV));

% create list of all faces 1, then 2, & 3
fnodes = [EToV(:,[1,2]);EToV(:,[2,3]);EToV(:,[3,1])];
fnodes = sort(fnodes,2)-1;

% set up default element to element and Element to faces connectivity
EToE= (1:K)'*ones(1,Nfaces);
EToF= ones(K,1)*(1:Nfaces);

% uniquely number each set of three faces by their node numbers 
id = fnodes(:,1)*Nnodes + fnodes(:,2)+1;
spNodeToNode=[id, (1:Nfaces*K)', EToE(:), EToF(:)];

% Now we sort by global face number.
sorted=sortrows(spNodeToNode,1);

% find matches in the sorted face list
[indices,dummy]=find( sorted(1:(end-1),1)==sorted(2:end,1) );

% make links reflexive 
matchL = [sorted(indices,:)   ;sorted(indices+1,:)];
matchR = [sorted(indices+1,:) ;sorted(indices,:)];

% insert matches
EToE(matchL(:,2)) = matchR(:,3);
EToF(matchL(:,2)) = matchR(:,4);

return;
