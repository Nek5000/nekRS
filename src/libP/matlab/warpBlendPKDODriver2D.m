clear  lebConstant
clear  condV


%% maximum polynomial degree
maxN = 12;

%% nodes used for Lebesgue calculation and plotting
plotN = 50;
[plotr,plots] = EquiNodes2D(plotN);
[plotr,plots] = xytors(plotr,plots);
plotNp = length(plotr);
plotTri = delaunay(plotr,plots);
figure(5)
for N=1:maxN

  %% number of modes and nodes
  Np = (N+1)*(N+2)/2;
  
  %% generate Warp & Blend Nodes for biunit right-angled triangle
  [r,s] = Nodes2D(N);
  [r,s] = xytors(r,s);

  %% build Vandermonde matrix using orthonormal PKDO basis
  V = Vandermonde2D(N, r, s);
  plotV = Vandermonde2D(N, plotr, plots);

  %% interpolation matrix
  Imatrix = plotV/V;

  %% compute condition number of V
  condV(N) = cond(V);

  %% sum up abs value of Lagrange basis functions at plot nodes
  lebFunction = sum(abs(Imatrix), 2);

  %% estimate of Lebesgue constant
  lebConstant(N) = max(lebFunction);

  %% plot lebFunction
  subplot(4,3,N); 
  trisurf(plotTri,plotr, plots, lebFunction); shading interp;
  axis([-1 1 -1 1 0 lebConstant(N)]); 
  xlabel('r');
  ylabel('s');
  [foo,ha] = suptitle('Lebesgue Function');
  set(ha,  'FontName', 'TimesNewRoman');
  set(gca, 'FontName', 'TimesNewRoman');
  pause(.05)
end
myprint('-dpng', 'Figures/warpBlendPKDOLebesgueFunction2D.png', '-portrait');

%% plot Lebesgue constant as a function of polynomial degree
figure(2)
clf; plot(1:maxN, lebConstant, 'k-r*'); axis([1 maxN 1 max(lebConstant)])
ha = title('Lebesgue Constant'); xlabel('Polynomial Degree N');
set(gca, 'FontName', 'TimesNewRoman'); set(ha, 'FontName', 'TimesNewRoman');
myprint('-dpdf', 'Figures/warpBlendPKDOLebesgueConstant2D.pdf', '-portrait')

%% plot Vandermonde condition number as a function of polynomial degree
figure(3)
clf; plot(1:maxN, condV, 'k-r*'); axis([1 maxN 1 max(condV)])
ha = title('cond(Vandermonde Matrix)'); xlabel('Polynomial Degree N');
set(gca, 'FontName', 'TimesNewRoman'); set(ha, 'FontName', 'TimesNewRoman');
myprinfigure(2)


