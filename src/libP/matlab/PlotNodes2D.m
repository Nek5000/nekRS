
function PlotNodes2D(N, col)

% add lines between nodes
Nsegs = 100;
for i=1:N+1
 for j=1:Nsegs
   linesR(j, i) = -1 + (i-1)*2/N;
   linesS(j, i) = -1 + 0.5*(1-linesR(j,i))*(j-1)*2/(Nsegs-1);
 end
end

% map to warped element 
[linesX,linesY] = WarpBlendTransform2D(N, linesR, linesS);

% use rotations for lines in each direction
linesX1 = [linesX,...
cos(2*pi/3)*linesX - sin(2*pi/3)*linesY,...
cos(4*pi/3)*linesX - sin(4*pi/3)*linesY];

linesY1 = [linesY,...
sin(2*pi/3)*linesX + cos(2*pi/3)*linesY,...
sin(4*pi/3)*linesX + cos(4*pi/3)*linesY];

[linesX1,linesY1] = xytors(linesX1,linesY1);

% draw lines connecting nodes
clf
hold on;
for i=1:(N+1)*3
 plot( linesX1(:,i), linesY1(:,i),'k-', 'LineWidth', 1.75);
end

% plot spheres at X,Y
[sX,sY] = Nodes2D(N);
[sX,sY] = xytors(sX,sY);

% sphere surface coords (uses MATLAB intrinsic sphere function)
[xs,ys,zs] = sphere(20);
ra= 0.05;

% choose default sphere color
if(nargin==1) col = [.9 0 .9]; end

for n=1:length(sX)
	ha = surf(ra*xs+sX(n), ra*ys+sY(n), ra*zs); %, col*ones(size(xs))); %ones(size(xs))*col)
%shading interp
%  material metal; lighting gouraud; 
  
  set(ha, 'FaceColor', col);
  set(ha, 'EdgeColor', 'none')
  
%  lighting gouraud
  
end
hold off 

% format scene
axis off; view(2); axis equal;
% light; light
set(gcf, 'color', 'white')

fname = sprintf('warpBlendNodes2DN%02d.pdf', N);
myprint('-dpdf', fname)

syscmd = sprintf('pdfcrop %s', fname)
  system(syscmd)
