function [Dr,Ds,Dt] = DmatricesHex3D(N,r,s,t,V)

%% function [Dr,Ds,Dt] = DmatricesHex3D(N,r,s,t,V)
%% Purpose : Initialize the (r,s,t) differentiation matrices
%%	    on the simplex, evaluated at (r,s,t) at order N

[Vr, Vs, Vt] = GradVandermondeHex3D(N, r, s, t);
Dr = Vr/V; Ds = Vs/V; Dt = Vt/V;

ids = find(abs(Dr)<1e-13); Dr(ids) = 0;
Dr = sparse(Dr);

ids = find(abs(Ds)<1e-13); Ds(ids) = 0;
Ds = sparse(Ds);

ids = find(abs(Dt)<1e-13); Dt(ids) = 0;
Dt = sparse(Dt);

return;
