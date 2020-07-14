function [P] = Simplex2DP(a,b,i,j);

% function [P] = Simplex2DP(a,b,i,j);
% Purpose : Evaluate 2D orthonormal polynomial
%           on simplex at (a,b) of order (i,j).

h1 = JacobiP(a,0,0,i); h2 = JacobiP(b,2*i+1,0,j);
P = sqrt(2.0)*h1.*h2.*(1-b).^i;
return;
