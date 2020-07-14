function [EToV] = delaunayOriented2D(VX, VY)

EToV = delaunay(VX, VY);
K = size(EToV, 1);
N = 1;
[r,s] = Nodes2D(N);
[r,s] = xytors(r,s);
V = Vandermonde2D(N, r, s);
[Dr, Ds] = Dmatrices2D(N, r, s, V);

for loop=1:3

  va = EToV(:,1)'; vb = EToV(:,2)'; vc = EToV(:,3)';
  x = 0.5*(-(r+s)*VX(va)+(1+r)*VX(vb)+(1+s)*VX(vc));
  y = 0.5*(-(r+s)*VY(va)+(1+r)*VY(vb)+(1+s)*VY(vc));
  
  %% build C0 stiffness
  [rx,sx,ry,sy,J] = GeometricFactors2D(x,y,Dr,Ds);
  
  %% remove slithers and reorient negative elements
  cnt = 1;
  for k=1:K
    if(min(abs(J(:,k)))>1e-6)
      if(min(J(:,k))<0)
	EToV(cnt,[1 3 2]) = EToV(k,:);
      else
	EToV(cnt,[1 2 3]) = EToV(k,:);
      end
      cnt = cnt + 1;
    end
  end
  K = cnt-1;
  EToV = EToV(1:K,:);

end

