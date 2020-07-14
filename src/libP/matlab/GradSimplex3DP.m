function [V3Dr, V3Ds, V3Dt] = GradSimplex3DP(a,b,c,id,jd,kd)

% function [V3Dr, V3Ds, V3Dt] = GradSimplex3DP(a,b,c,id,jd,kd)
% Purpose: Return the derivatives of the modal basis (id,jd,kd) 
%          on the 3D simplex at (a,b,c)

fa = JacobiP(a,0,0,id);           dfa = GradJacobiP(a,0,0,id);
gb = JacobiP(b,2*id+1,0,jd);      dgb = GradJacobiP(b,2*id+1,0,jd);
hc = JacobiP(c,2*(id+jd)+2,0,kd); dhc = GradJacobiP(c,2*(id+jd)+2,0,kd);

% r-derivative
V3Dr = dfa.*(gb.*hc);
if(id>0);    V3Dr = V3Dr.*((0.5*(1-b)).^(id-1));    end
if(id+jd>0); V3Dr = V3Dr.*((0.5*(1-c)).^(id+jd-1)); end

% s-derivative 
V3Ds = 0.5*(1+a).*V3Dr;
tmp = dgb.*((0.5*(1-b)).^id);
if(id>0);   tmp = tmp+(-0.5*id)*(gb.*(0.5*(1-b)).^(id-1));  end
if(id+jd>0) tmp = tmp.*((0.5*(1-c)).^(id+jd-1));            end
tmp = fa.*(tmp.*hc);
V3Ds = V3Ds+tmp;

% t-derivative 
V3Dt = 0.5*(1+a).*V3Dr+0.5*(1+b).*tmp;
tmp = dhc.*((0.5*(1-c)).^(id+jd));
if(id+jd>0)
  tmp = tmp-0.5*(id+jd)*(hc.*((0.5*(1-c)).^(id+jd-1))); 
end
tmp = fa.*(gb.*tmp); tmp = tmp.*((0.5*(1-b)).^id);
V3Dt = V3Dt+tmp;

% normalize
V3Dr = V3Dr*(2^(2*id+jd+1.5));
V3Ds = V3Ds*(2^(2*id+jd+1.5));
V3Dt = V3Dt*(2^(2*id+jd+1.5));
return
