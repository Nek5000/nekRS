% V = VDM.  Vrst = deriv matrices.  V1-4 barycentric derivs.  ids =
% permutation to match equispaced node ordering on tet. 

function [V Vr Vs Vt V1 V2 V3 V4 id] = bern_basis_tet(N,r,s,t)

%[L1 L2 L3 L4] = rsttobary(r,s,t);

L1 = -(1+r+s+t)/2; L2 = (1+r)/2; L3 = (1+s)/2; L4 = (1+t)/2;

dL1r = -.5; dL2r = .5; dL3r = 0; dL4r = 0;
dL1s = -.5; dL2s = 0; dL3s = .5; dL4s = 0;
dL1t = -.5; dL2t = 0; dL3t = 0; dL4t = .5;

sk = 1;
for l = 0:N
    for k = 0:N-l
        for j = 0:N-k-l;
            i = N-j-k-l;
            C=factorial(N)/(factorial(i)*factorial(j)*factorial(k)*factorial(l));
%             C = 1;
            V(:,sk) = C*(L1.^i).*(L2.^j).*(L3.^k).*(L4.^l);
            
            dL1 = C*i*(L1.^(i-1)).*(L2.^j).*(L3.^k).*(L4.^l);
            dL2 = C*j*(L1.^(i)).*(L2.^(j-1)).*(L3.^k).*(L4.^l);
            dL3 = C*k*(L1.^(i)).*(L2.^j).*(L3.^(k-1)).*(L4.^l);
            dL4 = C*l*(L1.^(i)).*(L2.^j).*(L3.^k).*(L4.^(l-1));
            if i==0
                dL1 = zeros(size(dL1));
            end
            if j==0
                dL2 = zeros(size(dL1));
            end
            if k ==0
                dL3 = zeros(size(dL1));
            end
            if l ==0
                dL4 = zeros(size(dL1));
            end
            Vr(:,sk) = dL1.*dL1r + dL2.*dL2r + dL3.*dL3r + dL4.*dL4r;
            Vs(:,sk) = dL1.*dL1s + dL2.*dL2s + dL3.*dL3s + dL4.*dL4s;
            Vt(:,sk) = dL1.*dL1t + dL2.*dL2t + dL3.*dL3t + dL4.*dL4t;
            V1(:,sk) = dL1;
            V2(:,sk) = dL2;
            V3(:,sk) = dL3;
            V4(:,sk) = dL4;
            sk = sk + 1;
        end
    end
end

function bi = bern(N,i,r)

bi = nchoosek(N,i)*(r.^i).*(1-r).^(N-i);

function dbi = d_bern(N,i,r)

dbi = nchoosek(N,i)*r.^(i - 1).*(1 - r).^(N - i - 1).*(i - N*r);

