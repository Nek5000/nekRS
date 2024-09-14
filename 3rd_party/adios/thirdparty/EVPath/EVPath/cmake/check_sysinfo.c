#include <sys/sysinfo.h>
int main() 
{
    struct sysinfo z;
    long sz = sysinfo(&z);
    unsigned int t = z.mem_unit;
    return (sz == 0) ? 0 : 1;
}
