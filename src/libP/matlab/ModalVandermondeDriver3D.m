tol = 1e-12;

for N=1:10

[r,s,t] = Nodes3D(N);
[r,s,t] = xyztorst(r,s,t);
Np = length(r);

[vertV,edgeV,intfaceV, intV] = ModalABCVandermonde3D(N, r, s, t);
mV = [vertV,edgeV,intfaceV,intV];
%mV = ModalVandermonde3D(N, r, s, t);

V = Vandermonde3D(N, r, s, t);
[Dr,Ds,Dt] = Dmatrices3D(N, r, s, t, V);
MM = inv(V*transpose(V));

Srr = Dr'*MM*Dr;
Srs = Dr'*MM*Ds + Ds'*MM*Dr;
Srt = Dr'*MM*Dt + Dt'*MM*Dr;

Sss = Ds'*MM*Ds;
Sst = Ds'*MM*Dt + Dt'*MM*Ds;

Stt = Dt'*MM*Dt;

mSrr = mV'*Srr*mV;
mSrs = mV'*Srs*mV;
mSrt = mV'*Srt*mV;
mSss = mV'*Sss*mV;
mSst = mV'*Sst*mV;
mStt = mV'*Stt*mV;

nnz(abs(mSrr)>tol)/(Np*Np)
nnz(abs(mSrs)>tol)/(Np*Np)
nnz(abs(mSrt)>tol)/(Np*Np)
nnz(abs(mSss)>tol)/(Np*Np)
nnz(abs(mSst)>tol)/(Np*Np)
nnz(abs(mStt)>tol)/(Np*Np)


subplot(3,3,1); spy(abs(mSrr)>tol)
subplot(3,3,2); spy(abs(mSrs)>tol)
subplot(3,3,3); spy(abs(mSrt)>tol)
subplot(3,3,5); spy(abs(mSss)>tol)
subplot(3,3,6); spy(abs(mSst)>tol)
subplot(3,3,9); spy(abs(mStt)>tol)

8*(nnz(abs(mSrr)>tol) + nnz(abs(mSrs)>tol) + nnz(abs(mSrt)>tol) + nnz(abs(mSss)>tol) + nnz(abs(mSst)>tol) + nnz(abs(mStt)>tol))/1024

[N,Np*Np*6*8/1024, Np*Np*8/1024]



end
