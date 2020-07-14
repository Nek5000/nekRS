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

NLambda = 2; maxLambda  = sqrt(1/2);

allv = [];
alld = [];
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
  [v,d] = eig(OP);
  d = diag(d);
  alld = [alld,d];
  allv = [allv,v];
subplot(1,3,1) 
hold on
  ha = plot(real(d), imag(d), 'r.', 'markersize', 10);
  set(ha, 'color', [Lambda/maxLambda,0,1-Lambda/maxLambda])
  xlabel('Re(lambda)', 'interpreter', 'tex')
  ylabel('Im(lambda)', 'interpreter', 'tex')
  axis([-40 .5 -200 200])
hold off
  drawnow
  pause(.05)

end

 hold off

 while(1)
   subplot(1,3,1)
   [gx,gy] =  ginput(1);
   
   [foo,ids] = min(abs(alld(:)-gx-i*gy));
   
   mode = reshape(allv(:,ids), Np, Nelements, 2);
   subplot(1,3,2)
   plot(x, real(mode(:,:,1)), 'r-', x, imag(mode(:,:,1)), 'b-')
   subplot(1,3,3)
   plot(x, real(mode(:,:,2)), 'r-', x, imag(mode(:,:,2)), 'b-')
   
				%fname = sprintf('scanAcousticsDG1DSpectrumN%d.pdf', N);
				%myprint('-dpdfwrite', fname);
				%system(sprintf('pdfcrop %s', fname))
 end
  
