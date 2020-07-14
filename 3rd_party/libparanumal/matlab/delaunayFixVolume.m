function tets = delaunayFixVolume(r,s,t)

tets = delaunay(r,s,t);

K = size(tets,1);

voltets = [];

sk = 1;
tol = 1e-10;
for k=1:K
    id1 = tets(k,1);  ids = tets(k,2:end);
    edges = [r(ids)-r(id1),s(ids)-s(id1),t(ids)-t(id1)];
    vol = abs(det(edges)/6);
    if(vol>tol)
        tetids(sk) = k;
        voltets = [voltets;vol];
        sk = sk + 1;
    end
end
vol = sum(voltets);
tets = tets(tetids,:);

