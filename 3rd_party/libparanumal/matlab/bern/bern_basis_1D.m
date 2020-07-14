function [V Vr] = bern_basis_1D(N,r)

r = (1+r)/2; % convert to unit
for i = 0:N
    V(:,i+1) = bern_1D(N,i,r);
    Vr(:,i+1) = d_bern(N,i,r)*.5; % change of vars
end

function bi = bern_1D(N,i,r)

bi = nchoosek(N,i)*(r.^i).*(1-r).^(N-i);

function dbi = d_bern(N,i,r)

if (i==0)
    dbi = -N*(1-r).^(N-1);
elseif (i==N)
    dbi = N*(r.^(N-1));
else
    dbi = nchoosek(N,i)*r.^(i - 1).*(1 - r).^(N - i - 1).*(i - N*r);
end

