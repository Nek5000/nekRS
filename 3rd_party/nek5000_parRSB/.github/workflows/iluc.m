load 'A.txt'
load 'LL.txt'
load 'UU.txt'

A = spconvert(A);
LL = spconvert(LL);
UU = spconvert(UU);
err = norm(A - LL * UU, Inf);
printf('LU error: %f\n', err);
assert(err < 1e-8);
