%% limits for domain [a,b]
a = -2; b = 3;
%% number of elements
Nelements = 3;
%% mesh vertices
VX = linspace(a, b, Nelements+1);
%% polynomial degree & number of nodes per element
N = 6;
Np = N+1;
%% load GLL interpolant nodes
r = JacobiGL(0,0,N);
%% Vandermonde
V = Vandermonde1D(N, r);
%% physical coordinates of all interpolation nodes
x = 0.5*(1-r)*VX(1:end-1) + 0.5*(1+r)*VX(2:end);
%% function with some random jitters
f = sin(pi*x) + rand(Np, Nelements);

%% plotting stuff
Nplot = 100;
rplot = transpose(linspace(-1,1, Nplot+1));
Vplot = Vandermonde1D(N, rplot);
xplot = 0.5*(1-rplot)*VX(1:end-1) + 0.5*(1+rplot)*VX(2:end);

%% GLL to plot nodes interpolation matrix 
Iplot = Vplot/V;

%% plot interpolant
figure(1)
plot(xplot, Iplot*f)
hold on
plot(x, f, 'k.',  'markersize', 15)
hold off
axis([-2,3, -1, 2])

%% plot reconstruction
figure(2)
fhat2 = f;
fhat2(:,1) = 0.5*(f(end,1)+f(1,2));
fhat2(:,3) = 0.5*(f(1,3)+f(end,2));
plot(xplot, Iplot*fhat2)
hold on
plot(x(Np:end-Np+1), f(Np:end-Np+1), 'k.',  'markersize', 15)
hold off
axis([-2,3, -1, 2])

figure(1)
xlabel('x')
myprint('-dpdfwrite', 'discontinuousInterpolantGLL1D.pdf')

figure(2)
xlabel('x')
myprint('-dpdfwrite', 'discontinuousInterpolantReconstructionGLL1D.pdf')

%% number of face nodes and faces
Nfp = 1;
Nfaces = 2;

% reference element collocation differentiation matrix
D = Dmatrix1D(N, r, V);

%% volume geometric factors
xr = D*x;
rx = 1./xr;
J = xr;

%% surface geometric factors
sJ = ones(Nfp*Nfaces, Nelements); %% surface Jacobian
nx = [-1;1]*ones(1,Nelements);  %% unit normal

% Lagrange coefficients of end nodes
L0 = zeros(Np,1);
L1 = zeros(Np,1);
L0(1)  = 1.0; 
L1(Np) = 1.0;

% reference element lift matrix
LIFT = V*(V'*[L0,L1]);


%% distributional derivative: volume terms
dfdr = D*f;
locdfdx = rx.*dfdr;

%% indices of end points
cnt = 1;
vmapM = zeros(Nfp*Nfaces, Nelements);
vmapP = zeros(Nfp*Nfaces, Nelements);
for e=1:Nelements
  vmapM(cnt) = 1  + (e-1)*Np;
  if(cnt>1)
   vmapP(cnt) = Np + (e-2)*Np;
  else
   vmapP(1) = 1;
  end
  cnt = cnt+1;
  vmapM(cnt) = Np + (e-1)*Np;
  if(cnt<Nelements*Nfp*Nfaces)
   vmapP(cnt) = 1 +    (e)*Np;
  else
   vmapP(cnt) = cnt;
  end
  cnt = cnt+1;
end
  
%% jumps in f
df = zeros(Nfp*Nfaces,Nelements);
fstar = 0.5*(f(vmapP)+f(vmapM));
df(:) = fstar-f(vmapM);

%% add lift of jump to elemental derivatives

dfdx = locdfdx + (1./J).*(LIFT*(sJ.*nx.*df));

figure(3)
dfdxplot = Iplot*dfdx;
locdfdxplot = Iplot*locdfdx;
plot(xplot, Iplot*dfdx, 'r-', 'linewidth', 2)
hold on
plot(xplot, Iplot*locdfdx, 'b-', 'linewidth', 2);
plot([x(1);x(end)],[0,0], 'k--')
hold off

xlabel('x')

title('Distributional derivative')
myprint('-dpdfwrite', 'distributionalDerivativeAverageGLL1D.pdf')
