function [V2Dr,V2Ds] = GradVandermonde2D(N,r,s)

% function [V2Dr,V2Ds] = GradVandermonde2D(N,r,s)
% Purpose : Initialize the gradient of the modal basis (i,j)
%		at (r,s) at order N

V2Dr = zeros(length(r),(N+1)*(N+2)/2);
V2Ds = zeros(length(r),(N+1)*(N+2)/2);

% find tensor-product coordinates
[a,b] = rstoab(r,s);

% Initialize matrices
sk = 1;
for i=0:N
  for j=0:N-i
    [V2Dr(:,sk),V2Ds(:,sk)] = GradSimplex2DP(a,b,i,j);
    sk = sk+1;
  end
end
return;
