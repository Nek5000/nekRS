
function [urhs,prhs] = acousticsRHS1D(A, Dr, drdx, LIFT, vmapM, vmapP, mapB, nx, Lambda, u, p)
  
  %% compute traces 
  uefM = u(vmapM);
  uefP = u(vmapP);
  pefM = p(vmapM);
  pefP = p(vmapP);

  pefP(mapB) = -pefM(mapB);
  
  %% compute intermediate state on element traces
  uefS = 0.5*(uefP+uefM) + (Lambda^2)*nx.*(A(1,1)*(uefP-uefM) + A(1,2)*(pefP-pefM));
  pefS = 0.5*(pefP+pefM) + (Lambda^2)*nx.*(A(2,1)*(uefP-uefM) + A(2,2)*(pefP-pefM));

  %% right hand sides
  dudr = Dr*u;
  dpdr = Dr*p;

  urhs = drdx.*(A(1,1)*dudr + A(1,2)*dpdr + ...
		LIFT*(nx.*(A(1,1)*(uefS-uefM) + A(1,2)*(pefS-pefM))));

  prhs = drdx.*(A(2,1)*dudr + A(2,2)*dpdr + ...
		LIFT*(nx.*(A(2,1)*(uefS-uefM) + A(2,2)*(pefS-pefM))));
