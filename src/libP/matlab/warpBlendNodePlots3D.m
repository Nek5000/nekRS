
%% maximum polynomial degree
maxN = 2;

%% nodes used for Lebesgue calculation and plotting
plotN = 12;
[plotr,plots,plott] = EquiNodes3D(plotN);
[plotr,plots,plott] = xyztorst(plotr,plots,plott);

Nplot = 10;
[plotr,plots,plott] = meshgrid(linspace(-1,1,Nplot));

plotNp = length(plotr);
plotTri = delaunay3(plotr(:),plots(:),plott(:));
figure(1)
for N=1:maxN

  %% number of modes and nodes
  Np = (N+1)*(N+2)*(N+3)/2;
  
  %% generate Warp & Blend Nodes for biunit right-angled triangle
  [r,s,t] = Nodes3D(N);
  [r,s,t] = xyztorst(r,s,t);

  %% build Vandermonde matrix using orthonormal PKDO basis
  V = Vandermonde3D(N, r, s, t);
  plotV = Vandermonde3D(N, plotr(:), plots(:), plott(:));

  %% interpolation matrix
  Imatrix = plotV/V;
  clf

  nidx = 1;
  for k=0:N
    for j=0:N-k
      for i=0:N-j-k
      
	%% plot lebFunction
        offsetr = 4*i;
        offsets = 4*j;
        offsett = 4*k;
	Im = reshape(Imatrix(:,nidx), size(plotr));

        cnt = 1;
	for c=linspace(-1,1,Nplot)
	  for b=linspace(-1,1,Nplot)
	    for a=linspace(-1,1,Nplot)
	      if( 1+a+b+c>0)
%  	      Im(cnt) = 0;
              end
              cnt = cnt+1;
	    end
	  end
	end

	iso = linspace(0.1, 1, 9);
	hold on
	for n=1:length(iso)
	ha = isosurface(plotr, plots, plott,Im,iso(n)); 
	end
	hold off
	drawnow;
	pause(.05)
	box on
	grid off
	xlabel('r');
	ylabel('s');
	zlabel('t');
%	tname = sprintf('Degree %d Lagrange Basis Functions', N);
%	[foo,ha] = suptitle(tname);
%	set(ha,  'FontName', 'TimesNewRoman');
%        set(gca, 'FontName', 'TimesNewRoman');
        drawnow;
        pause(.05)
        nidx = nidx+1;
       end
     end
  end
  
  fname = sprintf('Figures/warpBlendNodeN%02d.pdf', N);
  myprint('-dpdfwrite', fname);
 hold off;		  
end
