#include <stdio.h>
#include <ctype.h>
#include "cblas.h"
#include "cblas_f77.h"

#define XerblaStrLen 6
#define XerblaStrLen1 7

void
#ifdef HAS_ATTRIBUTE_WEAK_SUPPORT
__attribute__((weak))
#endif
F77_xerbla_base
#ifdef F77_CHAR
(F77_CHAR F77_srname, void *vinfo
#else
(char *srname, void *vinfo
#endif
#ifdef BLAS_FORTRAN_STRLEN_END
, size_t len
#endif
)
{
#ifdef F77_CHAR
   char *srname;
#endif

   char rout[] = {'c','b','l','a','s','_','\0','\0','\0','\0','\0','\0','\0'};

   int *info=vinfo;
   int i;

   extern int CBLAS_CallFromC;

#ifdef F77_CHAR
   srname = F2C_STR(F77_srname, XerblaStrLen);
#endif

   if (CBLAS_CallFromC)
   {
      for(i=0; i != XerblaStrLen; i++) rout[i+6] = tolower(srname[i]);
      rout[XerblaStrLen+6] = '\0';
      cblas_xerbla(*info+1,rout,"");
   }
   else
   {
      fprintf(stderr, "Parameter %d to routine %s was incorrect\n",
              *info, srname);
   }
}
