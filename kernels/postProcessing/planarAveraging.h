int mod1(int i, int n) {
  if(!i) return 0;
  return (i+n-1)%n + 1;
}

void get_exyz(int* ex, int* ey, int* ez, int eg, int nelx, int nely, int nelz)
{
  // convert 0-based indexing to 1-based indexing
  eg++;

  // compute based on 1-based indexing
  *ex = mod1(eg, nelx);
  *ey = 1 + (mod1(eg, nelx*nely) - 1)/nelx;
  *ez = 1 + (eg-1)/(nelx*nely);

  // shift back to 0-based indexing
  (*ex)--;
  (*ey)--;
  (*ez)--;

}