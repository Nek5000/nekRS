function [VrQ, VsQ] = VandermondeQuad2D(NQ, rQ, sQ)

  sk = 1;
  VrQ = zeros(length(rQ), (NQ+1)*(NQ+1));
  VsQ = zeros(length(rQ), (NQ+1)*(NQ+1));
  for i=0:NQ
    for j=0:NQ
      VrQ(:,sk) = GradJacobiP(rQ, 0, 0, i).*JacobiP(sQ, 0, 0, j);
      VsQ(:,sk) = JacobiP(rQ, 0, 0, i).*GradJacobiP(sQ, 0, 0, j);
      sk = sk+1;
    end
  end
  
