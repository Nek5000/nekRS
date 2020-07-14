function plotLiftPolynomials2DV2(N)

Globals2D;

Np = (N+1)*(N+2)/2;
Nfp = N+1;
[r,s] = Nodes2D(N);
[r,s] = xytors(r,s);
faceNodes = [find(abs(s+1)<1e-7),find(abs(r+s)<1e-7),find(abs(1+r)<1e-7)];

V = Vandermonde2D(N, r, s);

LIFT=Lift2D(N,faceNodes,r,s);

[rplot,splot] = EquiNodes2D(25);
[rplot,splot] = xytors(rplot,splot);

Vplot = Vandermonde2D(N, rplot, splot);
Iplot = Vplot/V;
triplot = delaunay(rplot,splot);
LIFTplot = Iplot*LIFT;

Nwin  = Nfp;
ids = reshape(1:Nwin*Nwin, Nwin, Nwin)';

clf
for n=1:Nfp
subplot(3,Nwin, n)
trisurf(triplot, rplot, splot, LIFTplot(:,n));
%trisurf(triplot, rplot, splot, log10(abs(LIFTplot(:,n)./Iplot(:,Fmask(n,1)))))

shading interp
view(3)
axis off
axis square
box on
end

for n=1:Nfp
subplot(3, Nwin, n+Nwin);
trisurf(triplot, rplot, splot, LIFTplot(:,1*Nfp+n))
shading interp
view(3)
axis off
axis square
box on

end


for n=1:Nfp
subplot(3,Nwin, n+2*Nwin);
trisurf(triplot, rplot, splot, LIFTplot(:,2*Nfp+n))
shading interp
view(3)
axis off
axis square
box on
end


fname = sprintf('trianglesLiftPolynomials2DN%d.pdf', N);

myprint('-dpdfwrite', fname)
system(sprintf('pdfcrop %s', fname))
