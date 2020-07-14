function [cV,cMM,cSrr,cSrs,cSss,stackedNz] = GenModalOps(N) 

  [r,s] = Nodes2D(N);
  [r,s] = xytors(r,s);
  Np = length(r);
  
  V = Vandermonde2D(N, r, s);
  [Dr, Ds] = Dmatrices2D(N, r, s, V);
  M = inv(V*transpose(V));
  invV = inv(V);

  %% continuous (sparse) basis
  cV = ModalVandermonde2D(N, r, s);
  
  %% change normalization
  cMM = cV'*M*cV;

  %cDMM = diag(diag(cMM));
  %cV = cV/sqrt(cDMM)
  %cMM = cV'*M*cV;

  %% compute derivative matrices
  cVr = Dr*cV;
  cVs = Ds*cV;

  %% build matrices
  cSrr = transpose(cVr)*M*cVr;

  cSrs = transpose(cVr)*M*cVs;
  cSrs = cSrs+transpose(cSrs);

  cSss = transpose(cVs)*M*cVs;

  %zero out very small entries
  cTol = 1e-10;
  cV(abs(cV)<cTol) = 0;
  cMM(abs(cMM)<cTol) = 0;
  cSrr(abs(cSrr)<cTol) = 0;
  cSrs(abs(cSrs)<cTol) = 0;
  cSss(abs(cSss)<cTol) = 0;
  
  %% build mask [ 1 when any of these three have a non-zero entry ]
  cMask = (abs(cSrr)>cTol) | (abs(cSrs)>cTol) | (abs(cSss)>cTol);

  maxNzPerRow = 0;		     
  for r=1:Np
   maxNzPerRow = max(maxNzPerRow, length(find(cMask(r,:))));		     
  end
		     
  stackedNz = zeros(Np, maxNzPerRow);
  for r=1:Np
   rowNz = sort(find(cMask(r,:)));		     
   rowNz = [rowNz,zeros(1,maxNzPerRow-length(rowNz))];
   stackedNz(r,:) = rowNz;
  end
	
end