function [V1D] = Vandermonde1D(N,xp)

% function [V1D] = Vandermonde1D(N,xp)
% Purpose : Initialize the 1D Vandermonde Matrix.
%	    V_{ij} = phi_j(xp_i);

V1D = zeros(length(xp),N+1);
for j=1:N+1
    V1D(:,j) = JacobiP(xp, 0, 0, j-1);
end;
return
