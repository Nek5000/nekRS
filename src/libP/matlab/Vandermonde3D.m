function [V3D] = Vandermonde3D(N, r, s, t);

% function [V3D] = Vandermonde3D(N, r, s, t);
% Purpose : Initialize the 3D Vandermonde Matrix, V_{ij} = phi_j(r_i, s_i, t_i);

V3D = zeros(length(r),(N+1)*(N+2)*(N+3)/6);

% Transfer to (a,b) coordinates
[a, b, c] = rsttoabc(r, s, t);

% build the Vandermonde matrix
sk = 1;

for i=0:N % old ordering
  for j=0:N - i
    for k=0:N - i - j
      V3D(:,sk) = Simplex3DP(a,b,c,i,j,k);
      sk = sk+1;
    end
  end
end
return;
