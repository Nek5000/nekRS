#include <sys/types.h>
#include <sys/sysctl.h>
main()
{
  long val = 0;
  static int mib[2] = {CTL_HW, HW_CPU_FREQ};
  size_t vlen = sizeof(long);
  sysctl(mib,2,&val,&vlen,NULL,0);
}
