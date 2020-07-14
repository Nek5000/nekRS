
%% maximum polynomial degree
maxN = 6;

%% nodes used for Lebesgue calculation and plotting
plotN = 12;
[plotr,plots] = EquiNodes2D(plotN);
[plotr,plots] = xytors(plotr,plots);
plotNp = length(plotr);
plotTri = delaunay(plotr,plots);
figure(1)
for N=1:maxN

  %% number of modes and nodes
  Np = (N+1)*(N+2)/2;
  
  %% generate Warp & Blend Nodes for biunit right-angled triangle
  [r,s] = Nodes2D(N);
  [r,s] = xytors(r,s);

  %% build Vandermonde matrix using orthonormal PKDO basis
  V = Vandermonde2D(N, r, s);
  plotV = Vandermonde2D(N, plotr, plots);

  %% interpolation matrix
  Imatrix = plotV/V;
  
  clf
  cnt = 1+(N+1)*N
  nidx = 1;
  for j=0:N
    for i=0:N-j
      node = cnt+i;
      %% plot lebFunction
      subplot(N+1,N+1,node);
      trisurf(plotTri,plotr, plots, Imatrix(:,nidx)); shading interp;
      drawnow;
      pause(.05)
      axis([-1 1 -1 1 -1.5 1.5], 'nolabel')
      box on
      grid off
      xlabel('r');
      ylabel('s');
      tname = sprintf('Degree %d Lagrange Basis Functions', N);
      [foo,ha] = suptitle(tname);
      set(ha,  'FontName', 'TimesNewRoman');
      set(gca, 'FontName', 'TimesNewRoman');
      drawnow;
      pause(.05)
      nidx = nidx+1;
    end
    cnt = cnt-(N+1)
  end
  
  fname = sprintf('Figures/warpBlendNodeN%02d.pdf', N);
  myprint('-dpdfwrite', fname);
		  
end
