%% maximum polynomial degree
maxN = 15;

%% number of plotting nodes [ also used for estimating Lebesgue constant ]
Nplot = 1000;

%% dimension of node set to test
Ntest = 10000;

%% polynomial degree test
figure(1)
for N=1:maxN
  %% number of points for unique interpolant
  Np = N+1

  %% build chebyshev  nodes with monomial basis
  rtest = chebyshevNodes1D(Ntest);
  rplot = transpose(linspace(-1,1,Nplot));
  
  %% build monomial Vandermonde
  Vtest = legendrePolynomials1D(N, rtest);
  Vplot = legendrePolynomials1D(N, rplot);

  %% use rank revealing qr to select "quasi-optimally orthogonal nodes"
  [Q,R,P] = qr(Vtest, '0');

  %% choose nodes using the qr pivot
  r = rtest(P(1:Np));
  V = Vtest(:,P(1:Np));

  %% interpolation matrix for nodes=>plot nodes
  Imatrix = Vtest/V;

  %% compute condition number of V
  condV(N) = cond(V);

  %% sum up abs value of Lagrange basis functions at plot nodes
  lebFunction = sum(abs(Imatrix), 2);

  %% estimate of Lebesgue constant
  lebConstant(N) = max(lebFunction);

  %% plot lebFunction
  subplot(3,5,N); plot(rtest,lebFunction, 'r-'); 
  axis([-1 1 0 lebConstant(N)]); xlabel('r');
  [foo,ha] = suptitle('Lebesgue Function');
  set(ha,  'FontName', 'TimesNewRoman');
  set(gca, 'FontName', 'TimesNewRoman');
end
myprint('-dpdf', 'Figures/rrqrLegendreLebesgueFunction1D.pdf', '-landscape');

%% plot Lebesgue constant as a function of polynomial degree
figure(2)
clf; plot(1:maxN, lebConstant, 'k-r*'); axis([1 maxN 1 max(lebConstant)])
ha = title('Lebesgue Constant'); xlabel('Polynomial Degree N');
set(gca, 'FontName', 'TimesNewRoman'); set(ha, 'FontName', 'TimesNewRoman');
myprint('-dpdf', 'Figures/rrqrLegendreLebesgueConstant1D.pdf', '-landscape')

%% plot Vandermonde condition number as a function of polynomial degree
figure(3)
clf; plot(1:maxN, condV, 'k-r*'); axis([1 maxN 1 max(condV)])
ha = title('cond(Vandermonde Matrix)'); xlabel('Polynomial Degree N');
set(gca, 'FontName', 'TimesNewRoman'); set(ha, 'FontName', 'TimesNewRoman');
myprint('-dpdf', 'Figures/rrqrLegendreVandermondeConditionNumber1D.pdf', '-landscape')

%% plot Basis functions
figure(4)
clf; plot(rplot, Vplot); 
ha = title('Basis functions(N=15)'); xlabel('r');
set(gca, 'FontName', 'TimesNewRoman'); set(ha, 'FontName', 'TimesNewRoman');
myprint('-dpdf', 'Figures/rrqrLegendreBasisFunctions1D.pdf', '-landscape')
