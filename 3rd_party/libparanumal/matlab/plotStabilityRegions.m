%% LOW STORAGE EXPLICIT RUNGE KUTTA


%% LSERK3
rk3a = [0, -5/9, -153/128];
rk3b = [1/3, 15/16, 8/15];


%% LSERK4
rk4a = [            0.0 ...
        -567301805773.0/1357537059087.0 ...
        -2404267990393.0/2016746695238.0 ...
        -3550918686646.0/2091501179385.0  ...
        -1275806237668.0/842570457699.0];
rk4b = [ 1432997174477.0/9575080441755.0 ...
         5161836677717.0/13612068292357.0 ...
         1720146321549.0/2090206949498.0  ...
         3134564353537.0/4481467310338.0  ...
         2277821191437.0/14882151754819.0];

[x,y] = meshgrid(-5:.01:5,-5:.01:5);
alpha = x + i*y;

z3 = ones(size(x));
resz3 = z3;
for INTRK3 = 1:3
  resz3 = rk3a(INTRK3)*resz3+ alpha.*z3;
  z3 = z3+rk3b(INTRK3)*resz3;
end

z4 = ones(size(x));
resz4 = z4;
for INTRK4 = 1:5    
  % 
  resz4 = rk4a(INTRK4)*resz4+ alpha.*z4;
  z4 = z4+rk4b(INTRK4)*resz4;
end

%% plot level set |z3|==1
%% plot level set |z3|==1
C3 = contourc(x,y,abs(z3),[0.;1.]); %, 'b-', 'linewidth', 2)
C4 = contourc(x,y,abs(z4),[0.;1.]); %, 'r-', 'linewidth', 2)

plot(C3(1,2:end)/3,C3(2,2:end)/3, 'r-', 'linewidth', 2)
hold on
plot(C4(1,2:end)/5,C4(2,2:end)/5, 'b-', 'linewidth', 2)
hold off
axis equal

%% ADAMS BASHFORTH STABILITY REGIONS
form{1} = 'r--'
form{2} = 'k--'
form{3} = 'b--'
form{4} = 'r--'

hold on;
for tstep=3:4
  
switch(tstep)
  case 1
    a = 1;
    a(1) = 1;   
  case 2
    a = [3 -1]/2;
  case 3
    a = [23 -16 5]/12;
  case {4,5}
    a = [17 9 1 -7]/20;
  case 6
    a = [1925 572 -316 -604 -337 485]/1680;
end

t = linspace(-pi, pi, 1000);

et = exp(i*t);

N = length(a);

deno = zeros(size(et));
for n=1:N
  deno = deno + a(n)*(et.^(N-n));
end

alpha = (et.^N - et.^(N-1))./deno;

if(tstep==4)
  ids = find(real(alpha)>-1.25);
alpha = alpha(ids);
end
plot(real(alpha), imag(alpha), form{tstep}, 'linewidth', 2)
%fill(real(alpha), imag(alpha),[tstep,tstep,tstep]/10)
end
axis equal


%% all legends
ha = legend('LSERK3', 'LSERK4', 'AB3', 'AB4');
set(ha, 'location', 'north')

axis equal

axis([-1.25 .25 -.8 .8])

xlabel('Re(lambda)')
ylabel('Im(lambda)')
title('Absolute stability regions')

myprint('-dpdfwrite', 'ABvRKStabiityRegions.pdf')
system('pdfcrop ABvRKStabiityRegions.pdf')



