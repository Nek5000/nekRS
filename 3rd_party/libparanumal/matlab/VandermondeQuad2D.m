function [V2D] = VandermondeQuad2D(N, r, s);

% function [V2D] = VandermondeQuad2D(N, r, s);
% Purpose : Initialize the 2D Vandermonde Matrix.
%           V_{ij} = phi_j(r_i, s_i);


V2D = zeros(length(r),(N+1)*(N+2)/2);

% build the Vandermonde matrix
sk = 1;
for i=0:N
  for j=0:N
    V2D(:,sk) = JacobiP(r, 0, 0, i).*JacobiP(s, 0, 0, j);
    sk = sk+1;
  end
end
return;
