
function writeFloatMatrix(file, A, label)
  
  fprintf(file, '******************************************\n');
  fprintf(file, '%s\n', label);
  fprintf(file, '%d %d\n', size(A,1), size(A,2));
  for r=1:size(A,1)
    for c=1:size(A,2)
      fprintf(file, '%17.15f ', A(r,c));
    end
    fprintf(file, '\n');
  end
  
  
