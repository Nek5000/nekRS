form{1} = 'r-'
form{2} = 'k-'
form{3} = 'b-'
form{4} = 'g--'

hold on;
for tstep=1:4
  
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
axis([-2 0.25 -1 1])
%axis([-2.5 .5 -1.5 1.5])

ha = legend('AB1', 'AB2', 'AB3', 'AB4'); %, "location", "northeastoutside");
xlabel('Re(lambda)')
ylabel('Im(lambda)')


myprint('-dpdfwrite', 'adamsBashforthStabiityRegions.pdf')
