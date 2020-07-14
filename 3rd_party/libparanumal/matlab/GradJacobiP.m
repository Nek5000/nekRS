function [dP] = GradJacobiP(z, alpha, beta, N);

% function [dP] = GradJacobiP(z, alpha, beta, N);
% Purpose: Evaluate the derivative of the orthonormal Jacobi
%	   polynomial of type (alpha,beta)>-1, at points x
%          for order N and returns dP[1:length(xp))]

dP = zeros(length(z), 1);
if(N == 0)
  dP(:,:) = 0.0; 
else
  dP = sqrt(N*(N+alpha+beta+1))*JacobiP(z,alpha+1,beta+1, N-1);
end;
return
