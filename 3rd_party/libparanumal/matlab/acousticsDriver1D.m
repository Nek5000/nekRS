%% interpolant degree
N = 6;

%% domain limits
xmin = -1;
xmax = +1;

%% number of elements
Nelements = 10;

%% build mesh and node info
[Np, EToV, VX, V, Dr, x, drdx, invM, vmapM, vmapP, mapB, nx, LIFT] = ...
    startUp1D(N, Nelements, xmin, xmax);

%% system
A = [ 0 1; 1 0];

%% Penalty parameter
Lambda = sqrt(.5); 

%% Euler forward time stepping
cfl = .2;
dt = cfl*min(min(1./drdx))/( (N+1)^2);
finalTime = 5;
Nsteps = ceil(finalTime/dt);
dt = finalTime/Nsteps;

%% initial solution
%u = cos(pi*x);
u = exp(-20*x.^2);
p = zeros(Np,Nelements);

for tstep=1:Nsteps

  %% compute RHS 
  [urhs,prhs] = acousticsRHS1D(A, Dr, drdx, LIFT, vmapM, vmapP, mapB, nx, Lambda, u, p);

  %% update u and p using Euler forward
  u = u + dt*urhs;
  p = p + dt*prhs;

  %% plot solution
  if(~mod(tstep, 40)||tstep==Nsteps)
    time = dt*tstep;

    clf
    plot(x, u,'r-');
    hold on;
    plot(x, p, 'b-'); 
    hold off;
    axis([-1 1 -1 1])
    title(sprintf('time = %f (u red, p blue)', time))

    drawnow; pause(.025);
		  
    maxError = max(max(abs(u - cos(pi*x)*cos(pi*time))))
  end
end


