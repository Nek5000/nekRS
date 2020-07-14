
function r = chebyshevNodes1D(N)

  %% Chebyshev node distribution on the bi-unit interval
  r = cos(0.5*pi*(2*(0:N)+1)/(N+1))';
