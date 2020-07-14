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

%% build DG operator (column by column) for a range of Lambda values
OP = zeros(2*Np*Nelements);
q = zeros(Np,Nelements*2);

NLambda = 100; maxLambda  = sqrt(16);
clf;
hold on
for Lambda=linspace(0, maxLambda, NLambda)

  %% for each node for each element for each field
  for n=1:2*Np*Nelements
    %% set one node value for one field to 1
    q(n) = 1;
    u = q(:,1:Nelements);
    p = q(:,Nelements+1:2*Nelements);

    %% operate on the vector with one non-zero in
    %% the n'th location to obtain the n'th column
    [urhs,prhs] = ...
	acousticsRHS1D(A, Dr, drdx, LIFT, vmapM, vmapP, mapB, nx, Lambda, u, p);
    OP(:,n) = [urhs(:);prhs(:)];
    q(n) = 0;
  end
  %% plot the spectrum for this Lambda
  d = eig(OP);
  ha = plot(real(d), imag(d), 'r.', 'markersize', 10);
  set(ha, 'color', [Lambda/maxLambda,0,1-Lambda/maxLambda])
  xlabel('Re(lambda)', 'interpreter', 'tex')
  ylabel('Im(lambda)', 'interpreter', 'tex')
%  axis([-5 .5 -200 200])
  drawnow
  pause(.05)

end

 hold off

axis equal

fname = sprintf('scanAcousticsDG1DSpectrumN%d.pdf', N);
myprint('-dpdfwrite', fname);
system(sprintf('pdfcrop %s', fname))

  
