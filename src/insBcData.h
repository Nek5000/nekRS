// used in boundary device functions
struct bcData
{
   int idM;
   int fieldOffset;
   int id;

   int scalarId;

   dfloat time;
   dfloat x, y, z;
   dfloat nx, ny, nz;

   dfloat uM, vM, wM;
   dfloat uP, vP, wP;
   dfloat uxP, uyP, uzP;
   dfloat vxP, vyP, vzP;
   dfloat wxP, wyP, wzP;

   dfloat pM;
   dfloat pP;

   @globalPtr const dfloat* wrk;

   dfloat sM, sP, sF; 
};
