function [Dr,Ds,Dt] = Dmatrices3D(N,r,s,t,V)

% function [Dr,Ds,Dt] = Dmatrices3D(N,r,s,t,V)
% Purpose : Initialize the (r,s,t) differentiation matrices
%	        on the simplex, evaluated at (r,s,t) at order N

[Vr, Vs, Vt] = GradVandermonde3D(N, r, s, t);
Dr = Vr/V; Ds = Vs/V; Dt = Vt/V;
return;
