
#include "cthread.h"
#include "io_malloc.h"
/* 
 * This uses a bit of a trick to gain convenience, so I'll explain.  The IO
 * routines need to allocate memory.  Sometimes they might need to use Cthread
 * memory allocation, but most times using malloc() is OK.  We arrange things
 * so that which is used can be chosen at link time.  This works because of
 * the way the linker operates.  To use Cthreads memory allocation, the user
 * must call IO_use_Cthread_malloc() before any other IOfile calls.  The call
 * to this routine causes the linker to include the object file from this C
 * file in the executable.  In this case, the program must be linked with
 * cthreads.  Alternatively, to use malloc(), the users *does not* include a
 * reference to the routine io_use_Cthread_malloc() and *does not* have to
 * link with Cthreads.  In this case, this object file is not included by the
 * linker and there should be no "memory_alloc undefined symbol" errors from
 * the linker.
 */

void *
io_Cthread_malloc(size)
int size;
{
    RESULT result;
    void *tmp;

    if ((result = memory_alloc((memory_t *) & tmp, size,
			       N_CURRENT)) != T_SUCCEED) {
	cthread_perror("memory allocation failed", result);
    }
    return tmp;
}

void *
io_Cthread_realloc(ptr, size)
void *ptr;
int size;
{
    RESULT result;

    if ((result = memory_realloc((memory_t *) & ptr, size)) != T_SUCCEED) {
	cthread_perror("memory allocation failed", result);
    }
    return ptr;
}

void
io_Cthread_free(ptr)
void *ptr;
{
    memory_free((memory_t) ptr);
}

extern
void
IO_use_Cthread_malloc()
{
    io_malloc_routine = io_Cthread_malloc;
    io_realloc_routine = io_Cthread_realloc;
    io_free_routine = io_Cthread_free;
}
