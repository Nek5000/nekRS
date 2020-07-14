%% maximum polynomial degree
minN = 1;
maxN = 15;

%% number of plotting nodes [ also used for estimating Lebesgue constant ]
Nplot = 1000;

%% dimension of node set to test
Ntest = 8000;

%% polynomial degree test
figure(1)
for N=minN:maxN
  %% number of points for unique interpolant
  Np = N+1

  %% build chebyshev  nodes with monomial basis
  rplot = chebyshevNodes1D(Nplot);
  Vplot = legendrePolynomials1D(N, rplot);

  %% scan over blending factor
  lebesgueConstant = 1e9;
  optimalAlpha = -1;
  for alpha=linspace(0, .5, Ntest)
    ralpha = blend1D(N, alpha);
    Valpha = legendrePolynomials1D(N, ralpha);
    Ialpha = Vplot/Valpha;
    lebesgueConstantAlpha = max(sum(abs(Ialpha),2));
    if(lebesgueConstantAlpha<lebesgueConstant)
      lebesgueConstant = lebesgueConstantAlpha;
      optimalAlpha = alpha;
      saveAlpha(N) = optimalAlpha;
    end
  end

  r = blend1D(N, optimalAlpha);
  V = legendrePolynomials1D(N, r);

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
myprint('-dpdf', 'Figures/blendLegendreLebesgueFunction1D.pdf', '-landscape');

%% plot Lebesgue constant as a function of polynomial degree
figure(2)
clf; plot(1:maxN, lebConstant, 'k-r*'); axis([1 maxN 1 max(lebConstant)])
ha = title('Lebesgue Constant'); xlabel('Polynomial Degree N');
set(gca, 'FontName', 'TimesNewRoman'); set(ha, 'FontName', 'TimesNewRoman');
myprint('-dpdf', 'Figures/blendLegendreLebesgueConstant1D.pdf', '-landscape')

%% plot Vandermonde condition number as a function of polynomial degree
figure(3)
clf; plot(1:maxN, condV, 'k-r*'); axis([1 maxN 1 max(condV)])
ha = title('cond(Vandermonde Matrix)'); xlabel('Polynomial Degree N');
set(gca, 'FontName', 'TimesNewRoman'); set(ha, 'FontName', 'TimesNewRoman');
myprint('-dpdf', 'Figures/blendLegendreVandermondeConditionNumber1D.pdf', '-landscape')

%% plot Basis functions
figure(4)
clf; plot(rplot, Vplot); 
ha = title('Basis functions(N=15)'); xlabel('r');
set(gca, 'FontName', 'TimesNewRoman'); set(ha, 'FontName', 'TimesNewRoman');
myprint('-dpdf', 'Figures/blendLegendreBasisFunctions1D.pdf', '-landscape')

%% myprint out Lebesgue constants
lebConstant
