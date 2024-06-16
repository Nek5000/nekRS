void timeInterpWeights(dfloat xx, const dfloat *x, int n, dfloat *c)
{
  dfloat c1 = 1.0;
  dfloat c4 = x[0] - xx;

  for (int j = 0; j <= n; j++) {
    c[j] = 0.0;
  }
  c[0] = 1.0;

  for (int i = 1; i <= n; i++) {
    dfloat c2 = 1.0;
    dfloat c5 = c4;
    c4 = x[i] - xx;
    for (int j = 0; j <= i - 1; j++) {
      dfloat c3 = x[i] - x[j];
      c2 = c2 * c3;
      c[i] = -c1 * c5 * c[i - 1] / c2;
      c[j] = c4 * c[j] / c3;
    }
    c1 = c2;
  }
}