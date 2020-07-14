function [X,Y,Z] = EquiNodes3D(N)

% function [X,Y,Z] = EquinNdes3D(N)
% Purpose: compute the equidistributed nodes on the reference tetrahedron

% total number of nodes
Np = (N+1)*(N+2)*(N+3)/6;

% 2) create equidistributed nodes on equilateral triangle
X = zeros(Np,1); Y = zeros(Np,1); Z = zeros(Np,1); 

sk = 1;
for n=1:N+1
  for m=1:N+2-n
    for q=1:N+3-n-m
      X(sk) = -1 + (q-1)*2/N; Y(sk) = -1 + (m-1)*2/N; Z(sk) = -1 + (n-1)*2/N;
      sk = sk+1;
    end
  end
end
return;
