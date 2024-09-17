#include <unistd.h>
/*May have return value backwards.  if returns error, exits with 1.*/
int main() 
{
  long val;
  val = (long) sysconf(_SC_NPROCESSORS_ONLN);
  return ( val == -1 ) ? 1 : 0;
}
