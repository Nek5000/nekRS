function [a,b] = rstoab(r,s)

% function [a,b] = rstoab(r,s)
% Purpose : Transfer from (r,s) -> (a,b) coordinates in triangle

Np = length(r); a = zeros(Np,1);
for n=1:Np
  if(s(n) ~= 1)
     a(n) = 2*(1+r(n))/(1-s(n))-1;
  else a(n) = -1; end
end
b = s;
return;
