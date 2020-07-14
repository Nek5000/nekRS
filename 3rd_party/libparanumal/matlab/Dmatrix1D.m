function [Dr] = Dmatrix1D(N,r,V)

% function [Dr] = Dmatrix1D(N,r,V)
% Purpose : Initialize the (r) differentiation matrices
%	    on the interval, evaluated at (r) at order N

Vr = GradVandermonde1D(N, r);
Dr = Vr/V;
return
