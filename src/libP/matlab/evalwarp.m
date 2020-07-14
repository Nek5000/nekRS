function warp = evalwarp(p, xnodes, xout)

% function warp = evalwarp(p, xnodes, xout)
% Purpose: compute one-dimensional edge warping function

warp = zeros(size(xout));

for i=1:p+1
  xeq(i) = -1 + 2*(p+1-i)/p;
end

for i=1:p+1
  d = (xnodes(i)-xeq(i));
  for j=2:p
    if(i~=j)
	d = d.*(xout-xeq(j))/(xeq(i)-xeq(j));
    end
  end
  
  if(i~=1)
    d = -d/(xeq(i)-xeq(1));
  end

  if(i~=(p+1))
    d = d/(xeq(i)-xeq(p+1));
  end

  warp = warp+d;
end
return;
