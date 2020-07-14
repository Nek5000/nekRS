
%% return Legendre polynomials up to degree N
function [D,V] = legendreDifferentiationMatrix1D(N, r)

  %% evaluate Vandermonde and gradient of Vandermonde
  %% using Legendre polynomial basis
  Vr = zeros(length(r), N+1);
  V = zeros(length(r), N+1);

  V(:,1) = 1;
  Vr(:,1) = 0;

  V(:,2) = r;
  Vr(:,2) = 1;

  for id=3:N+1
    m = id-2;
    %%  (m+1)*L_(m+1)(r)-(2m+1)r*L_m(r)+m*L_(m-1)(r)=0
    V(:,id) = (1/(m+1))*( (2*m+1)*r.*V(:,id-1) - m*V(:,id-2)); 
    %%  (m+1)*dL_(m+1)(r)-(2m+1)(L_m(r)+r*dL_m(r)) +m*dL_(m-1)(r)=0
    Vr(:,id) = (1/(m+1))*( (2*m+1)*(V(:,id-1) + r.*Vr(:,id-1)) - m*Vr(:,id-2)); 
  end
  
  %% interpolate then differentiate
  D = Vr/V;
