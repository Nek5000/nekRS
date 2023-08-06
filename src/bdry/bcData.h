// used in boundary device functions
struct bcData
{
  int idM;

  int fieldOffset;
  int id;

  dfloat time;
  dfloat x, y, z;
  dfloat nx, ny, nz;

  // tangential directions
  dfloat t1x, t1y, t1z;
  dfloat t2x, t2y, t2z;

  dfloat tr1, tr2;

  dfloat u, v, w;
  dfloat p;

  // interpolated velocity values
  dfloat uinterp, vinterp, winterp;

  int scalarId;
  dfloat s, flux;

  // interpolated scalar value
  dfloat sinterp;

  dfloat meshu, meshv, meshw;

  // properties
  dfloat trans, diff;

  @globalPtr const dfloat* usrwrk;
};
