#include <sys/time.h>
#include <sys/resource.h>
int main(int argc, char** argv)
{
    struct rlimit rlim;
    int ret = getrlimit(RLIMIT_MEMLOCK, &rlim);
    return (rlim.rlim_max < 102400);
}
