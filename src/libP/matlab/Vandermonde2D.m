function [V2D] = Vandermonde2D(N, r, s);

% function [V2D] = Vandermonde2D(N, r, s);
% Purpose : Initialize the 2D Vandermonde Matrix.
%           V_{ij} = phi_j(r_i, s_i);

V2D = zeros(length(r),(N+1)*(N+2)/2);

% Transfer to (a,b) coordinates
[a, b] = rstoab(r, s);

% build the Vandermonde matrix
sk = 1;
for i=0:N
  for j=0:N - i
    V2D(:,sk) = Simplex2DP(a,b,i,j);
    sk = sk+1;
  end
end
return;
