function A = myspconvert(A, nrows, ncols, tol)

  % function A = myspconvert(A, nrows, ncols, tol)
  % Purpose: convert a triplet sparse matrix into a Matlab sparse
  %          matrix and crop entries smaller than tol 

  ids = find(A(:,1)==0);
  if(~isempty(ids))
    A = A(1:min(ids)-1,:);
  end

  ids = find(abs(A(:,3))>tol);
  A = A(ids, :);

  % convert to sparse matrix
  A = spconvert(A);

  % pad to required size
  if(size(A,1)<nrows | size(A,2)<ncols)
    A(nrows,ncols) = 0;
  end
 
 
