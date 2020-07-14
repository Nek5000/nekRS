function [cubR,cubS,cubW, Ncub] = Cubature2D(Corder)

% function [cubR,cubS,cubW, Ncub] = Cubature2D(Corder)
% Purpose: provide multidimensional quadrature (i.e. cubature) 
%          rules to integrate up to Corder polynomials
     
CubatureData2D;

if(Corder<=28)
  cubR = cub2D{Corder}(:,1);
  cubS = cub2D{Corder}(:,2);
  cubW = cub2D{Corder}(:,3); 
else
  cubNA = ceil( (Corder+1)/2);
  [cubA,cubWA] = JacobiGQ(0,0, cubNA-1);
  cubNB = ceil( (Corder+1)/2);
  [cubB,cubWB] = JacobiGQ(1,0, cubNB-1);
  
  cubA = ones(cubNB,1)*cubA';
  cubB = cubB*ones(1,cubNA);

  cubR = 0.5*(1+cubA).*(1-cubB)-1;
  cubS = cubB;
  cubW = 0.5*cubWB*(cubWA');

  cubR = cubR(:);
  cubS = cubS(:);
  cubW = cubW(:);
end
Ncub = length(cubW(:));

return
