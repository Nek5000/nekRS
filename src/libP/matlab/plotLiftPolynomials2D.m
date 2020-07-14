Globals2D;
N =  8
Np = (N+1)*(N+2)/2;
Nfp = N+1;
[r,s] = Nodes2D(N);
[r,s] = xytors(r,s);
faceNodes = [find(abs(s+1)<1e-7),find(abs(r+s)<1e-7),find(abs(1+r)<1e-7)];

V = Vandermonde2D(N, r, s);

LIFT=Lift2D(N,faceNodes,r,s);

[rplot,splot] = EquiNodes2D(30);
[rplot,splot] = xytors(rplot,splot);

Vplot = Vandermonde2D(N, rplot, splot);
Iplot = Vplot/V;
triplot = delaunay(rplot,splot);
LIFTplot = Iplot*LIFT;

Nwin  = Nfp+2;
ids = reshape(1:Nwin*Nwin, Nwin, Nwin)';

clf
for n=1:Nfp
subplot(Nwin, Nwin, ids(Nwin,1+n)); 
trisurf(triplot, rplot, splot, LIFTplot(:,n))
shading interp
axis off
box on
view(2)
end

for n=1:Nfp
subplot(Nwin, Nwin, ids(Nwin-n,Nwin-n)); 
trisurf(triplot, rplot, splot, LIFTplot(:,1*Nfp+n))
shading interp
axis off
box on
view(2)
end


for n=1:Nfp
subplot(Nwin, Nwin, ids(Nwin-n,1)); 
trisurf(triplot, rplot, splot, LIFTplot(:,2*Nfp+n))
shading interp
axis off
box on
view(2)
end
