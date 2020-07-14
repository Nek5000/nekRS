function [Dr,Ds] = DmatricesQuad2D(N,r,s,V)

% function [Dr,Ds] = DmatricesQuad2D(N,r,s,V)
% Purpose : Initialize the (r,s) differentiation matrices
%	    on the simplex, evaluated at (r,s) at order N

[Vr, Vs] = GradVandermondeQuad2D(N, r, s);
Dr = Vr/V; Ds = Vs/V;

ids = find(abs(Dr)<1e-13); Dr(ids) = 0;
Dr = sparse(Dr);

ids = find(abs(Ds)<1e-13); Ds(ids) = 0;
Ds = sparse(Ds);

return;
