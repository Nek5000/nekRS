function [V Vr Vs V1 V2 V3 id] = bern_basis_tri(N,r,s)


[V Vr Vs V1 V2 V3] = bern_tri(N,r,s);

return
% 
% % use equivalence between W&B and equispaced nodes - get ordering
% [re se] = EquiNodes2D(N); [re se] = xytors(re,se);
% Ve = bern_tri(N,re,se);
% for i = 1:size(Ve,2)
%    [val iid] = max(Ve(:,i)); 
%    id(i) = iid;
% %    id(i) = i;
% end
% 
% [V Vr Vs V1 V2 V3] = bern_tri(N,r,s);
% V  = V(:,id);
% Vr = Vr(:,id);
% Vs = Vs(:,id);
% V1 = V1(:,id);
% V2 = V2(:,id);
% V3 = V3(:,id);

function [V Vr Vs VL1 VL2 VL3] = bern_tri(N,r,s)

% barycentric version
L1 = -(r+s)/2; L2 = (1+r)/2; L3 = (1+s)/2;
dL1r = -.5; dL2r = .5; dL3r = 0;
dL1s = -.5; dL2s = 0; dL3s = .5;

sk = 1;
for k = 0:N
    for j = 0:N-k
        i = N-j-k;
        C=factorial(N)/(factorial(i)*factorial(j)*factorial(k));
        V(:,sk) = C*(L1.^i).*(L2.^j).*(L3.^k);
        
        dL1 = C*i*(L1.^(i-1)).*(L2.^j).*(L3.^k);
        dL2 = C*j*(L1.^(i)).*(L2.^(j-1)).*(L3.^k);
        dL3 = C*k*(L1.^(i)).*(L2.^j).*(L3.^(k-1));
        if i==0
            dL1 = zeros(size(dL1));
        end
        if j==0
            dL2 = zeros(size(dL2));
        end
        if k == 0
            dL3 = zeros(size(dL3));
        end
        Vr(:,sk) = dL1.*dL1r + dL2.*dL2r + dL3.*dL3r;        
        Vs(:,sk) = dL1.*dL1s + dL2.*dL2s + dL3.*dL3s;
        
        VL1(:,sk) = dL1;
        VL2(:,sk) = dL2;
        VL3(:,sk) = dL3;
        sk = sk + 1;
    end
end

return



% [a b] = rstoab(r,s);
% a = .5*(a+1); % convert to unit
% b = .5*(b+1); % convert to unit

% be careful about derivatives at s = 1
a = .5*(r+1)./(1-.5*(s+1));
b = .5*(s+1);
a(abs(1-s)<1e-8) = 0;

dadr = 1./(1 - s);
dads = (r+1)./((1-s).^2);
dbds = .5;

sk = 1;
for i = 0:N
    for j = 0:N-i
        k = N-i-j;
        V(:,sk) = bern(N-k,i,a).*bern(N,k,b);
        Vr(:,sk) = d_bern(N-k,i,a).*bern(N,k,b).*dadr;
        Vs(:,sk) = d_bern(N-k,i,a).*bern(N,k,b).*dads + bern(N-k,i,a).*d_bern(N,k,b)*dbds;
        sk = sk + 1;
    end
end

function bi = bern(N,i,r)

bi = nchoosek(N,i)*(r.^i).*(1-r).^(N-i);

function dbi = d_bern(N,i,r)

dbi = nchoosek(N,i)*r.^(i - 1).*(1 - r).^(N - i - 1).*(i - N*r);

