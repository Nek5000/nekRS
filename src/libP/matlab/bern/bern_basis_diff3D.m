function [D0ids,D0vals,D1ids,D1vals,D2ids,D2vals,D3ids,D3vals] = bern_basis_diff3D(N,V,VB,V0,V1,V2,V3)

Np = (N+1)*(N+2)*(N+3)/6;

%entries are integers
D0 = VB\V0; D0(abs(D0)< 0.5) = 0; %all the entries are integers
D1 = VB\V1; D1(abs(D1)< 0.5) = 0;
D2 = VB\V2; D2(abs(D2)< 0.5) = 0;
D3 = VB\V3; D3(abs(D3)< 0.5) = 0;

%convert to sparse storage

%4 non-zeros per row
D0ids = zeros(Np,4);
D1ids = zeros(Np,4);
D2ids = zeros(Np,4);
D3ids = zeros(Np,4);

D0vals = zeros(Np,4);
D1vals = zeros(Np,4);
D2vals = zeros(Np,4);
D3vals = zeros(Np,4);

for i = 1:Np
    tmp = find(D0(i,:));
    D0vals(i,1:length(tmp)) = D0(i,tmp);
    tmp = tmp-1; % zero indexing
    if length(tmp) < 4
        tmp = [tmp zeros(1,4-length(tmp))];
    end
    D0ids(i,:) = tmp;

    tmp = find(D1(i,:));
    D1vals(i,1:length(tmp)) = D1(i,tmp);
    tmp = tmp-1; % zero indexing
    if length(tmp) < 4
        tmp = [tmp zeros(1,4-length(tmp))];
    end
    D1ids(i,:) = tmp;


    tmp = find(D2(i,:));
    D2vals(i,1:length(tmp)) = D2(i,tmp);
    tmp = tmp-1; % zero indexing
    if length(tmp) < 4
        tmp = [tmp zeros(1,4-length(tmp))];
    end
    D2ids(i,:) = tmp;
    
    tmp = find(D3(i,:));
    D3vals(i,1:length(tmp)) = D3(i,tmp);
    tmp = tmp-1; % zero indexing
    if length(tmp) < 4
        tmp = [tmp zeros(1,4-length(tmp))];
    end
    D3ids(i,:) = tmp;
end


return
