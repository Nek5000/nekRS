function [D0Tids,D0Tvals,D1Tids,D1Tvals,D2Tids,D2Tvals,D3Tids,D3Tvals] = bern_basis_transposediff3D(N,V,VB,V0,V1,V2,V3)

Np = (N+1)*(N+2)*(N+3)/6;

%entries are integers
D0T = (VB\V0)'; D0T(abs(D0T)< 0.5) = 0; %all the entries are integers
D1T = (VB\V1)'; D1T(abs(D1T)< 0.5) = 0;
D2T = (VB\V2)'; D2T(abs(D2T)< 0.5) = 0;
D3T = (VB\V3)'; D3T(abs(D3T)< 0.5) = 0;

%convert to sparse storage

%4 non-zeros per row
D0Tids = zeros(Np,4);
D1Tids = zeros(Np,4);
D2Tids = zeros(Np,4);
D3Tids = zeros(Np,4);

D0Tvals = zeros(Np,4);
D1Tvals = zeros(Np,4);
D2Tvals = zeros(Np,4);
D3Tvals = zeros(Np,4);

for i = 1:Np
    tmp = find(D0T(i,:));
    D0Tvals(i,1:length(tmp)) = D0T(i,tmp);
    tmp = tmp-1; % zero indexing
    if length(tmp) < 4
        tmp = [tmp zeros(1,4-length(tmp))];
    end
    D0Tids(i,:) = tmp;

    tmp = find(D1T(i,:));
    D1Tvals(i,1:length(tmp)) = D1T(i,tmp);
    tmp = tmp-1; % zero indexing
    if length(tmp) < 4
        tmp = [tmp zeros(1,4-length(tmp))];
    end
    D1Tids(i,:) = tmp;


    tmp = find(D2T(i,:));
    D2Tvals(i,1:length(tmp)) = D2T(i,tmp);
    tmp = tmp-1; % zero indexing
    if length(tmp) < 4
        tmp = [tmp zeros(1,4-length(tmp))];
    end
    D2Tids(i,:) = tmp;
    
    tmp = find(D3T(i,:));
    D3Tvals(i,1:length(tmp)) = D3T(i,tmp);
    tmp = tmp-1; % zero indexing
    if length(tmp) < 4
        tmp = [tmp zeros(1,4-length(tmp))];
    end
    D3Tids(i,:) = tmp;
end


return
