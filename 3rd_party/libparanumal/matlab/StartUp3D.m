% Purpose : Setup script, building operators, grid, metric,
%           and connectivity tables for 3D meshes of tetrahedra.

% Definition of constants
Np = (N+1)*(N+2)*(N+3)/6; Nfp = (N+1)*(N+2)/2; Nfaces=4; NODETOL = 1e-7;

% Compute nodal set
[x,y,z] = Nodes3D(N); [r,s,t] = xyztorst(x,y,z);

% Build reference element matrices
V = Vandermonde3D(N,r,s,t); invV = inv(V);
MassMatrix = invV'*invV;
[Dr,Ds,Dt] = Dmatrices3D(N, r, s, t, V);

% build coordinates of all the nodes
va = EToV(:,1)'; vb = EToV(:,2)'; vc = EToV(:,3)'; vd = EToV(:,4)';
x = 0.5*(-(1+r+s+t)*VX(va)+(1+r)*VX(vb)+(1+s)*VX(vc)+(1+t)*VX(vd));
y = 0.5*(-(1+r+s+t)*VY(va)+(1+r)*VY(vb)+(1+s)*VY(vc)+(1+t)*VY(vd));
z = 0.5*(-(1+r+s+t)*VZ(va)+(1+r)*VZ(vb)+(1+s)*VZ(vc)+(1+t)*VZ(vd));

% find all the nodes that lie on each edge
fmask1   = find( abs(1+t) < NODETOL)'; 
fmask2   = find( abs(1+s) < NODETOL)';
fmask3   = find( abs(1+r+s+t) < NODETOL)';
fmask4   = find( abs(1+r) < NODETOL)';
Fmask  = [fmask1;fmask2;fmask3;fmask4]';
Fx = x(Fmask(:), :); Fy = y(Fmask(:), :); Fz = z(Fmask(:), :);

% Create surface integral terms
LIFT = Lift3D(N, Fmask, r, s, t);

% calculate geometric factors
[rx,sx,tx,ry,sy,ty,rz,sz,tz,J] = GeometricFactors3D(x,y,z,Dr,Ds,Dt);

% calculate geometric factors
[nx, ny, nz, sJ] = Normals3D();
Fscale = sJ./(J(Fmask,:));

% Build connectivity matrix
[EToE, EToF] = tiConnect3D(EToV); 

% Build connectivity maps
[vmapM, vmapP, vmapB, mapB] = BuildMaps3D();

% Compute weak operators (could be done in preprocessing to save time)
[Vr, Vs, Vt] = GradVandermonde3D(N, r, s, t);
Drw = (V*Vr')/(V*V'); Dsw = (V*Vs')/(V*V'); Dtw = (V*Vt')/(V*V');
