
function [r,s,allIds] = GroupNodes2D(N)

  [r0,s0] = Nodes2D(N);
  Np = length(r0);
  
  theta = 2*pi/3;
  
  r1 = cos(theta)*r0 - sin(theta)*s0;
  s1 = sin(theta)*r0 + cos(theta)*s0;

  r2 = cos(2*theta)*r0 - sin(2*theta)*s0;
  s2 = sin(2*theta)*r0 + cos(2*theta)*s0;
  
  one = ones(Np,1);

  d1 = (r1*one' - one*r0').^2 + (s1*one' - one*s0').^2;
  d2 = (r2*one' - one*r0').^2 + (s2*one' - one*s0').^2;
  
  [mind1,ids1] = min(d1,[], 1);
  [mind2,ids2] = min(d2,[], 1);

  allIds = [1:Np; ids1; ids2];

  mask = zeros(Np,1);

  ids = [];
  for n=1:Np
    if(sum(mask(allIds(:,n)))==0)
      if(length(unique(allIds(:,n)))==1)
	newids  = allIds(1,n); %% just one new node
      else
	newids = allIds(:,n);
      end
      mask(newids) = 1;
      ids = [ids;newids];
    end
  end

  r = r0(ids);
  s = s0(ids);
  
  
