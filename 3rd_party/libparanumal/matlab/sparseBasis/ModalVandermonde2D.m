
function V = ModalVandermonde2D(N, r, s)

  Np = (N+1)*(N+2)/2;
  
  [a,b] = rstoab(r,s);

  cnt = 1;

  V = zeros(length(r), Np);

  V(:,cnt) = 0.25*(1-a).*(1-b); cnt = cnt + 1;
  V(:,cnt) = 0.25*(1+a).*(1-b); cnt = cnt + 1;
  V(:,cnt) = 0.5*(1+b);         cnt = cnt + 1;

  for i=0:N-2

    %% edge 1
    q = i+1;
    V(:,cnt) = 0.125*(1+a).*(1-b).*(1+b).*JacobiP(b, 1, 1, q-1); 
    cnt = cnt+1;

    %% edge 2
    q = i+1;
    V(:,cnt) = 0.125*(1-a).*(1-b).*(1+b).*JacobiP(b, 1, 1, q-1); 
    cnt = cnt+1;

    %% interior modes
    for j=1:i

      p = j;
      q = i+1-j;

      phia = 0.25*(1+a).*(1-a).*JacobiP(a, 1, 1, p-1);
      V(:,cnt) = 0.5*phia.*((0.5*(1-b)).^(p+1)).*(1+b).*JacobiP(b,2*p+1,1,q-1); 
      cnt = cnt+1;
    end

    %% edge 0
    p = i+1;
    V(:,cnt) = 0.25*(1+a).*(1-a).*JacobiP(a, 1, 1, p-1).*((0.5*(1-b)).^(p+1)); 
    cnt = cnt+1;

  end

  mismatch = cnt-Np-1;
  if(mismatch)
    mismatch
  end