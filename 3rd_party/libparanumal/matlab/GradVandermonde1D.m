function [DVr] = GradVandermonde1D(N,r)

% function [DVr] = GradVandermonde1D(N,r)
% Purpose : Initialize the gradient of the modal basis (i)
%		at (r) at order N

DVr = zeros(length(r),(N+1));

% Initialize matrix
for i=0:N
   [DVr(:,i+1)] = GradJacobiP(r,0,0,i);
end
return
