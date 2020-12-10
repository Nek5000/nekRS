// used in boundary device functions
struct bcData
{
  int idM;

  int fieldOffset;
  int id;

  dfloat time;
  dfloat x, y, z;
  dfloat nx, ny, nz;

  dfloat u, v, w;
  dfloat p;

  int scalarId;
  dfloat s, flux;

  @globalPtr const dfloat* wrk;
};
