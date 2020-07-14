%% maximum polynomial degree
maxN = 15;

%% number of plotting nodes [ also used for estimating Lebesgue constant ]
Nplot = 1000;

figure(1)
%% polynomial degree test
for N=1:maxN
  %% number of points for unique interpolant
  Np = N+1;

  %% build chebyshev  nodes with monomial basis
  r     = transpose(linspace(-1,1,Np));
  rplot = transpose(linspace(-1,1,Nplot));
  
  %% transform to Chebyshev
  r     = chebyshevNodes1D(N);
  rplot = chebyshevNodes1D(Nplot);

  %% build monomial Vandermonde
  V = zeros(Np, N+1);
  Vplot = zeros(Nplot+1, N+1);
  for n=0:N
    V(:,n+1) = r.^n;
    Vplot(:,n+1) = rplot.^n;
  end

  %% interpolation matrix for nodes=>plot nodes
  Imatrix = Vplot/V;

  %% compute condition number of V
  condV(N) = cond(V);

  %% sum up abs value of Lagrange basis functions at plot nodes
  lebFunction = sum(abs(Imatrix), 2);

  %% estimate of Lebesgue constant
  lebConstant(N) = max(lebFunction);

  %% plot lebFunction
  subplot(3,5,N); plot(rplot,lebFunction, 'r-'); 
  axis([-1 1 0 lebConstant(N)]); xlabel('r');
  [foo,ha] = suptitle('Lebesgue Function');
  set(ha,  'FontName', 'TimesNewRoman');
  set(gca, 'FontName', 'TimesNewRoman');
end
myprint('-dpdf', 'Figures/chebyshevMonomialsLebesgueFunction1D.pdf', '-landscape');

%% plot Lebesgue constant as a function of polynomial degree
figure(2)
clf; plot(1:N, lebConstant, 'k-r*'); axis([1 maxN 1 max(lebConstant)])
ha = title('Lebesgue Constant'); xlabel('Polynomial Degree N');
set(gca, 'FontName', 'TimesNewRoman'); set(ha, 'FontName', 'TimesNewRoman');
myprint('-dpdf', 'Figures/chebyshevMonomialsLebesgueConstant1D.pdf', '-landscape')

%% plot Vandermonde condition number as a function of polynomial degree
figure(3)
clf; semilogy(1:N, condV, 'k-r*'); axis([1 maxN 1 max(condV)])
ha = title('cond(Vandermonde Matrix)'); xlabel('Polynomial Degree N');
set(gca, 'FontName', 'TimesNewRoman'); set(ha, 'FontName', 'TimesNewRoman');
myprint('-dpdf', 'Figures/chebyshevMonomialsVandermondeConditionNumber1D.pdf', '-landscape')

%% plot Basis functions
figure(4)
clf; plot(rplot, Vplot); 
ha = title('Basis functions(N=15)'); xlabel('r');
set(gca, 'FontName', 'TimesNewRoman'); set(ha, 'FontName', 'TimesNewRoman');
myprint('-dpdf', 'Figures/chebyshevMonomialsBasisFunctions1D.pdf', '-landscape')
