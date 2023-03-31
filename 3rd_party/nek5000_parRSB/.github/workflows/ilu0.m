load 'A.txt';
load 'B.txt';
A = spconvert(A);
B = spconvert(B);

[L, U] = ilu(A);
n = size(A, 1);
I = speye(n);
err = norm(L + U - B - I, Inf);
printf('ILU err = %f', err);
assert(err < 1e-8);
