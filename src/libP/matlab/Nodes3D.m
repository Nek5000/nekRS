function [X,Y,Z] = Nodes3D(p)

% function [X,Y,Z] = Nodes3D(p)
% Purpose: compute Warp & Blend nodes
%  input:    p=polynomial order of interpolant
%  output: X,Y,Z vectors of node coordinates in equilateral tetrahedron

% choose optimized blending parameter
alphastore = [0;0;0;0.1002; 1.1332;1.5608;1.3413;1.2577;1.1603;...
	            1.10153;0.6080;0.4523;0.8856;0.8717;0.9655];
if(p<=15); alpha = alphastore(p) ; else;  alpha = 1. ; end 

% total number of nodes and tolerance
N = (p+1)*(p+2)*(p+3)/6; tol = 1e-10;

[r,s,t] = EquiNodes3D(p); % create equidistributed nodes 
L1 = (1+t)/2; L2 = (1+s)/2; L3 = -(1+r+s+t)/2; L4 = (1+r)/2;

% set vertices of tetrahedron
v1 = [-1, -1/sqrt(3), -1/sqrt(6)]; v2 = [ 1, -1/sqrt(3),-1/sqrt(6)]; 
v3 = [ 0,  2/sqrt(3), -1/sqrt(6)]; v4 = [ 0,  0,         3/sqrt(6)];

% orthogonal axis tangents on faces 1-4
t1(1,:) = v2-v1;          t1(2,:) = v2-v1;
t1(3,:) = v3-v2;          t1(4,:) = v3-v1;   
t2(1,:) = v3-0.5*(v1+v2); t2(2,:) = v4-0.5*(v1+v2);
t2(3,:) = v4-0.5*(v2+v3); t2(4,:) = v4-0.5*(v1+v3);  

for n=1:4 % normalize tangents 
   t1(n,:) = t1(n,:)/norm(t1(n,:)); t2(n,:) = t2(n,:)/norm(t2(n,:)); 
end  

% Warp and blend for each face (accumulated in shiftXYZ)
XYZ = L3*v1+L4*v2+L2*v3+L1*v4; % form undeformed coordinates
shift = zeros(size(XYZ));
for face=1:4 
  if(face==1); La = L1; Lb = L2; Lc = L3; Ld = L4; end;
  if(face==2); La = L2; Lb = L1; Lc = L3; Ld = L4; end;
  if(face==3); La = L3; Lb = L1; Lc = L4; Ld = L2; end; 
  if(face==4); La = L4; Lb = L1; Lc = L3; Ld = L2; end;

  % compute warp tangential to face
  [warp1 warp2] = WarpShiftFace3D(p, alpha, alpha, La, Lb, Lc, Ld); 
  
  blend = Lb.*Lc.*Ld;   % compute volume blending

  denom = (Lb+.5*La).*(Lc+.5*La).*(Ld+.5*La);   % modify linear blend
  ids = find(denom>tol);
  blend(ids) = (1+(alpha.*La(ids)).^2).*blend(ids)./denom(ids);
  
  % compute warp & blend
  shift = shift+(blend.*warp1)*t1(face,:) + (blend.*warp2)*t2(face,:);
  
  % fix face warp 
  ids = find(La<tol & ( (Lb>tol) + (Lc>tol) + (Ld>tol) < 3));
  shift(ids,:) = warp1(ids)*t1(face,:) + warp2(ids)*t2(face,:);
end
XYZ = XYZ + shift;
X = XYZ(:,1); Y = XYZ(:,2); Z = XYZ(:,3);
return;
