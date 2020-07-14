function [P] = Simplex3DP(a,b,c,i,j,k);

% function [P] = Simplex3DP(a,b,c,i,j,k);
% Purpose : Evaluate 3D orthonormal polynomial
%           on simplex at (a,b,c) of order (i,j,k).

h1 = JacobiP(a,0,0,i); h2 = JacobiP(b,2*i+1,0,j); h3 = JacobiP(c,2*(i+j)+2,0,k);
P = 2*sqrt(2)*h1.*h2.*((1-b).^i).*h3.*((1-c).^(i+j));
return;
