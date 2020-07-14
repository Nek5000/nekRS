function [r, s, t] = xyztorst(X, Y, Z)

% function [r,s,t] = xyztorst(x, y, z)
% Purpose : Transfer from (x,y,z) in equilateral tetrahedron
%           to (r,s,t) coordinates in standard tetrahedron

v1 = [-1,-1/sqrt(3), -1/sqrt(6)]; v2 = [ 1,-1/sqrt(3), -1/sqrt(6)];
v3 = [ 0, 2/sqrt(3), -1/sqrt(6)]; v4 = [ 0, 0/sqrt(3),  3/sqrt(6)];

% back out right tet nodes
rhs = [X';Y';Z'] - 0.5*(v2'+v3'+v4'-v1')*ones(1,length(X));
A = [0.5*(v2-v1)',0.5*(v3-v1)',0.5*(v4-v1)'];
RST = A\[rhs];
r = RST(1,:)'; s = RST(2,:)'; t = RST(3,:)';
return;
