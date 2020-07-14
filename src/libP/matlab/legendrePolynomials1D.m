
%% return Legendre polynomials up to degree N
function L = legendrePolynomials1D(N, r)

  L = zeros(length(r), N+1);
  L(:,1) = 1;
  L(:,2) = r;

  %%  (m+1)P_(m+1)(r)-(2m+1)rP_m(r)+mP_(m-1)(r)=0
  for id=3:N+1
    m = id-2;
    L(:,id) = (1/(m+1))*( (2*m+1)*r.*L(:,id-1) - m*L(:,id-2)); 
  end
  
