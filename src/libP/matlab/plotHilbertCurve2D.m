
function plotHilbertCurve2D(D, nums)

H = 0:2^D;
[I,J] = meshgrid(H);
H = IJtoH(D,I,J);

H = 0:2^(2*D)-1;
[I,J] = HtoIJ(D, H);
plot(I(:)+.5,J(:)+.5, 'b-');
hold on
plot3(I+.5,J+.5, .01+0*I, 'r.', 'markersize', 10);

if(nums=='y')
  fsize = 10
  if(D==4) fsize = 8 end
  
  for n=1:length(I(:))
    ha = text(I(n)+.3,J(n)+.2,sprintf('%03d',n-1), 'fontsize', fsize);
  end
end
hold off

axis equal
axis off

fname = sprintf('hilbertCurve2D%dB.pdf', D)
myprint('-dpdfwrite', fname)
scmd = sprintf('pdfcrop %s', fname);
system(scmd)

return 
end



%% adapted from https://en.wikipedia.org/wiki/Hilbert_curve
%% convert (I,J) to Hilbert index H

function H = IJtoH(D, I, J)
  H = zeros(size(I));
  
  for n=D-1:-1:1

    rI = bitand(I, 2^n) > 0;
    rJ = bitand(J, 2^n) > 0;

    H = H + 2^(2*n) * bitxor((3 * rI), rJ);
    [I,J] = rot(2^n, I, J, rI, rJ);
  end
  return 
end

%% convert Hilbert index H to (I,J)
function [I, J] = HtoIJ(D, H)

  I = zeros(size(H));
  J = zeros(size(H));
  for n=0:D-1
    rI = bitand(1, floor(H/2));
    rJ = bitand(1, bitxor(H,rI));
    [I,J] = rot(2^n, I, J, rI, rJ);
    I = I + 2^n * rI;
    J = J + 2^n * rJ;
    H = floor(H/4);
  end
  return 
end
   

%% rotate/flip a quadrant appropriately
function [I,J] = rot(n, I, J, rI, rJ)

  for m=1:length(I(:))

    if (rJ(m) == 0)
      if (rI(m) == 1) 
	I(m) = n-1 - I(m);
        J(m) = n-1 - J(m);
      end
      
      %% Swap x and y
      t = I(m);
      I(m) = J(m);
      J(m) = t;
      
    end
  end
  return 
end