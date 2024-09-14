#ifdef WITH_ADIOS2
#include <adios2.h>
#endif

#include "foo.h"

void foo(void)
{
#ifdef WITH_ADIOS2
    adios2::ADIOS ctx;
#endif
}
