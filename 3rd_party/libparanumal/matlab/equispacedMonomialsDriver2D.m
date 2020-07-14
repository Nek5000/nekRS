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
  
  %% generate equispaced Nodes for biunit right-angled triangle
  [r,s] = EquiNodes2D(N);
  [r,s] = xytors(r,s);

  %% build Vandermonde matrix using monomials
  V = zeros(Np, Np);
  plotV = zeros(plotNp, Np);
  
  cnt = 1;
  for m2=0:N
    for m1=0:N-m2
      V(:,cnt) = (r.^m1).*(s.^m2);
      plotV(:,cnt) = (plotr.^m1).*(plots.^m2);
      cnt = cnt+1;
    end
  end

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
myprint('-dpng', 'Figures/equispacedMonomialsLebesgueFunction2D.png', '-portrait');

%% plot Lebesgue constant as a function of polynomial degree
figure(2)
clf; semilogy(1:maxN, lebConstant, 'k-r*'); axis([1 maxN 1 max(lebConstant)])
ha = title('Lebesgue Constant'); xlabel('Polynomial Degree N');
set(gca, 'FontName', 'TimesNewRoman'); set(ha, 'FontName', 'TimesNewRoman');
myprint('-dpdf', 'Figures/equispacedMonomialsLebesgueConstant2D.pdf', '-portrait')

%% plot Vandermonde condition number as a function of polynomial degree
figure(3)
clf; semilogy(1:maxN, condV, 'k-r*'); axis([1 maxN 1 max(condV)])
ha = title('cond(Vandermonde Matrix)'); xlabel('Polynomial Degree N');
set(gca, 'FontName', 'TimesNewRoman'); set(ha, 'FontName', 'TimesNewRoman');
myprint('-dpdf', 'Figures/equispacedMonomialsVandermondeConditionNumber2D.pdf', '-portrait')
