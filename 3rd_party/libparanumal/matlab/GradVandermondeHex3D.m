function [Vr, Vs, Vt] = VandermondeHex3D(N, r, s, t)

  sk = 1;
  Vr = zeros(length(r), (N+1)*(N+1)*(N+1));
  Vs = zeros(length(r), (N+1)*(N+1)*(N+1));
  Vt = zeros(length(r), (N+1)*(N+1)*(N+1));
  for k=0:N
    for j=0:N
      for i=0:N
	Vr(:,sk) = GradJacobiP(r, 0, 0, i).*JacobiP(s, 0, 0, j).*JacobiP(t, 0, 0, k);
	Vs(:,sk) = JacobiP(r, 0, 0, i).*GradJacobiP(s, 0, 0, j).*JacobiP(t, 0, 0, k);
	Vt(:,sk) = JacobiP(r, 0, 0, i).*JacobiP(s, 0, 0, j).*GradJacobiP(t, 0, 0, k);
	sk = sk+1;
      end
    end
  end
  
