function [V3Dr,V3Ds,V3Dt] = GradVandermonde3D(N,r,s,t)

% function [V3Dr,V3Ds,V3Dt] = GradVandermonde3D(N,r,s,t)
% Purpose : Initialize the gradient of the modal basis (i,j,k) at (r,s,t) at order N

V3Dr = zeros(length(r),(N+1)*(N+2)*(N+3)/6);
V3Ds = zeros(length(r),(N+1)*(N+2)*(N+3)/6);
V3Dt = zeros(length(r),(N+1)*(N+2)*(N+3)/6);

% find tensor-product coordinates
[a,b,c] = rsttoabc(r,s,t);

% Initialize matrices

sk = 1;
for i=0:N
  for j=0:N-i
    for k=0:N-i-j
      [V3Dr(:,sk),V3Ds(:,sk),V3Dt(:,sk)] = GradSimplex3DP(a,b,c,i,j,k);
      sk = sk+1;
    end
  end
end
return;
