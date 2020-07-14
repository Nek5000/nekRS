function [V3D] = VandermondeHex2D(N, r, s, t);

% function [V3D] = VandermondeHex2D(N, r, s, t);
% Purpose : Initialize the 3D Vandermonde Matrix.
%           V_{ij} = phi_j(r_i, s_i, t_i);


V3D = zeros(length(r),(N+1)*(N+1)*(N+1));

% build the Vandermonde matrix
sk = 1;
for k=0:N
  for j=0:N
    for i=0:N
      V3D(:,sk) = JacobiP(r, 0, 0, i).*JacobiP(s, 0, 0, j).*JacobiP(t, 0, 0, k);
      sk = sk+1;
    end
  end
end
return;
