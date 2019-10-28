// used in boundary device functions
struct bcData
{
   int idM;
   int fieldOffset;
   int id;

   int voffset;
   int soffset;

   dfloat time;
   dfloat x, y, z;
   dfloat nx, ny, nz;

   dfloat uM, vM, wM;
   dfloat uP, vP, wP;
   dfloat uxM, uyM, uzM;
   dfloat vxM, vyM, vzM;
   dfloat wxM, wyM, wzM;
   dfloat uxP, uyP, uzP;
   dfloat vxP, vyP, vzP;
   dfloat wxP, wyP, wzP;

   dfloat pM;
   dfloat pP, pxP, pyP, pzP;

   @globalPtr dfloat* wrk;

   dfloat sM, sP, sF; 
};
