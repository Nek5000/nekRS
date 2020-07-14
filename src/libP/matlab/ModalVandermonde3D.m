function V = ModalVandermonde3D(N, r, s, t)

  Nmodes = (N+1)*(N+2)*(N+3)/6;

  V = zeros(length(r), Nmodes);

  sk = 1;

  % vertices 1-4
  V(:,sk) = (-1/2)*(1+r+s+t); sk = sk+1;
  V(:,sk) = ( 1/2)*(1+r);     sk = sk+1;
  V(:,sk) = ( 1/2)*(1+s);     sk = sk+1;
  V(:,sk) = ( 1/2)*(1+t);     sk = sk+1;

  % edge 1
  for i=1:N-1
    V(:,sk) = (1/4)*(1+r).*(-(1+r+s+t)).*OldJacobiP(r, 1, 1, i-1);
    sk = sk+1;
  end

  % edge 2
  for i=1:N-1
    V(:,sk) = (1/4)*(1+r).*(1+s).*OldJacobiP(s, 1, 1, i-1);
    sk = sk+1;
  end

  % edge 3
  for i=1:N-1
    V(:,sk) = (1/4)*(-(1+r+s+t)).*(1+s).*OldJacobiP(s, 1, 1, i-1);
    sk = sk+1;
  end
  
  % edge 4
  for i=1:N-1
    V(:,sk) = (1/4)*(-(1+r+s+t)).*(1+t).*OldJacobiP(t, 1, 1, i-1);
    sk = sk+1;
  end

  % edge 5
  for i=1:N-1
    V(:,sk) = (1/4)*(1+r).*(1+t).*OldJacobiP(t, 1, 1, i-1);
    sk = sk+1;
  end

  % edge 6
  for i=1:N-1
    V(:,sk) = (1/4)*(1+s).*(1+t).*OldJacobiP(t, 1, 1, i-1);
    sk = sk+1;
  end

  % face 1
  for i=1:N-1
    for j=1:N-1-i
      V(:,sk) = (1/8)*(1+r).*(1+s).*(1+r+s+t).*OldJacobiP(r, 1, 1, i-1).*OldJacobiP(s, 1, 1, j-1);
      sk = sk+1;
    end
  end

  % face 2
  for i=1:N-1
    for j=1:N-1-i
      V(:,sk) = (-1/8)*(1+r).*(1+t).*(1+r+s+t).*OldJacobiP(r, 1, 1, i-1).*OldJacobiP(t, 1, 1, j-1);
      sk = sk+1;
    end
  end


  % face 3
  for i=1:N-1
    for j=1:N-1-i
      V(:,sk) = (1/8)*(1+r).*(1+s).*(1+t).*OldJacobiP(s, 1, 1, i-1).*OldJacobiP(t, 1, 1, j-1);
      sk = sk+1;
    end
  end

  % face 4
  for i=1:N-1
    for j=1:N-1-i
      V(:,sk) = (-1/8)*(1+s).*(1+t).*(1+r+s+t).*OldJacobiP(r, 1, 1, i-1).*OldJacobiP(t, 1, 1, j-1);
      sk = sk+1;
    end
  end
      
  % interior
  for i=1:N-1
    for j=1:N-1-i
      for k=1:N-1-i-j
	[i-1+j-1+k-1]
	V(:,sk) = ...
	    (-1/16)*(1+r).*(1+s).*(1+t).*(1+r+s+t).*OldJacobiP(r, 1, 1, i-1).*OldJacobiP(s, 1, 1, j-1).*OldJacobiP(t, 1, 1, k-1);
	sk = sk+1
      end
    end
  end
  
