/*
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C test:
 *  Write a huge number of variables
 *  Then read them all and check if they are correct.
 *
 * How to run: mpirun -np <N> many_vars <nvars> <blocks per process> <steps>
 * Output: many_vars.bp
 *
 */
#include "adios2_c.h"
#include <errno.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef WIN32
#define strncasecmp _strnicmp
#endif

#ifdef DMALLOC
#include "dmalloc.h"
#endif

#define log(...)                                                                                   \
    fprintf(stderr, "[rank=%3.3d, line %d]: ", rank, __LINE__);                                    \
    fprintf(stderr, __VA_ARGS__);                                                                  \
    fflush(stderr);
#define printE(...)                                                                                \
    fprintf(stderr, "[rank=%3.3d, line %d]: ERROR: ", rank, __LINE__);                             \
    fprintf(stderr, __VA_ARGS__);                                                                  \
    fflush(stderr);

int NVARS = 1;
int NBLOCKS = 1;
int NSTEPS = 1;
int REDEFINE = 0; // 1: delete and redefine variable definitions at each step to
                  // test adios_delete_vardefs()
static const char FILENAME[] = "many_vars.bp";
#define VALUE(rank, step, block) (step * 10000 + 10 * rank + block)

/* Variables to write */
int *a2;

static const int ldim1 = 5;
static const int ldim2 = 5;
int gdim1, gdim2;
int offs1, offs2;

/* ADIOS2 variables to be accessed from multiple functions */
adios2_io *ioW, *ioR;
adios2_engine *engineW;
adios2_variable **varW; // array to hold all variable definitions

/* Variables to read */
int *r2;

MPI_Comm comm = MPI_COMM_WORLD;
int rank;
int size;
char **varnames;

void alloc_vars()
{
    int n, i;

    n = ldim1 * ldim2;
    a2 = (int *)malloc(n * sizeof(int));
    r2 = (int *)malloc(n * sizeof(int));
    varW = (adios2_variable **)malloc(NVARS * sizeof(adios2_variable *));

    varnames = (char **)malloc(NVARS * sizeof(char *));
    for (i = 0; i < NVARS; i++)
    {
        varnames[i] = (char *)malloc(16);
    }

    /* make varnames like v001,v002,.. */
    int digit = 1, d = 10;
    while (NVARS / d > 0)
    {
        d *= 10;
        digit++;
    }

    char fmt[16];
    snprintf(fmt, sizeof(fmt), "v%%%d.%dd", digit, digit);
    printf("fmt=[%s]\n", fmt);
    for (i = 0; i < NVARS; i++)
    {
        // sprintf(varnames[i], "v%-*d", digit, i);
        snprintf(varnames[i], 16, fmt, i);
    }
    printf("varname[0]=%s\n", varnames[0]);
    printf("varname[%d]=%s\n", NVARS - 1, varnames[NVARS - 1]);
}

void set_gdim()
{
    gdim1 = size * ldim1;
    gdim2 = NBLOCKS * ldim2;
}

void set_offsets(int block)
{
    offs1 = rank * ldim1;
    offs2 = block * ldim2;
}

void set_vars(int step, int block)
{
    int n, i;
    int v = VALUE(rank, step, block);

    set_offsets(block);

    n = ldim1 * ldim2;
    for (i = 0; i < n; i++)
        a2[i] = v;
}

void fini_vars()
{
    int i;
    free(a2);
    free(r2);
    free(varW);
    for (i = 0; i < NVARS; i++)
    {
        free(varnames[i]);
    }
    free(varnames);
}

void Usage()
{
    printf("Usage: many_vars <nvars> <nblocks> <nsteps> [redef]\n"
           "    <nvars>:   Number of variables to generate\n"
           "    <nblocks>: Number of blocks per process to write\n"
           "    <nsteps>:  Number of write cycles (to same file)\n"
           "    [redef]:   delete and redefine variables at every step\n");
}

void define_vars();
int write_file(int step);
int read_file();

int main(int argc, char **argv)
{
    int err, i;

    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (argc < 2)
    {
        Usage();
        return 1;
    }

    errno = 0;
    i = strtol(argv[1], NULL, 10);
    if (errno || i < 1)
    {
        printf("Invalid 1st argument %s\n", argv[1]);
        Usage();
        return 1;
    }
    NVARS = i;

    if (argc > 2)
    {
        i = strtol(argv[2], NULL, 10);
        if (errno || i < 1)
        {
            printf("Invalid 2nd argument %s\n", argv[2]);
            Usage();
            return 1;
        }
        NBLOCKS = i;
    }

    if (argc > 3)
    {
        i = strtol(argv[3], NULL, 10);
        if (errno || i < 1)
        {
            printf("Invalid 3rd argument %s\n", argv[3]);
            Usage();
            return 1;
        }
        NSTEPS = i;
    }

    if (argc > 4)
    {
        if (!strncasecmp(argv[4], "redef", 5))
        {
            printf("Delete and redefine variable definitions at each step.\n");
            REDEFINE = 1;
        }
    }

    alloc_vars();
    adios2_adios *adiosH = adios2_init_mpi(MPI_COMM_WORLD);
    ioW = adios2_declare_io(adiosH, "multiblockwrite"); // group for writing
    ioR = adios2_declare_io(adiosH, "multiblockread");  // group for reading
    set_gdim();

    engineW = adios2_open(ioW, FILENAME, adios2_mode_write);

    err = 0;
    for (i = 0; i < NSTEPS; i++)
    {
        if (!err)
        {
            if (i == 0 || REDEFINE)
            {
                printf("-- Define variables.\n");
                define_vars();
            }

            err = write_file(i);

            if (REDEFINE)
            {
                printf("-- Delete variable definitions.\n");
                adios2_remove_all_variables(ioW);
                adios2_remove_all_attributes(ioW);
            }
        }
    }
    adios2_close(engineW);

    if (!err)
        err = read_file();

    adios2_finalize(adiosH);
    fini_vars();
    MPI_Finalize();
    return err;
}

void define_vars()
{
    int i;

    size_t shape[2] = {gdim1, gdim2};
    size_t count[2] = {ldim1, ldim2};
    size_t start[2] = {0, 0};

    /* One variable definition for many blocks.
     * Offsets will change at writing for each block. */
    for (i = 0; i < NVARS; i++)
    {
        varW[i] = adios2_define_variable(ioW, varnames[i], adios2_type_int32_t, 2, shape, start,
                                         count, adios2_constant_dims_false);
    }
}

int write_file(int step)
{
    int block, v, i;
    double tb, te;
    size_t count[2] = {ldim1, ldim2};

    log("Write step %d to %s\n", step, FILENAME);
    tb = MPI_Wtime();

    adios2_step_status status;
    adios2_begin_step(engineW, adios2_step_mode_append, -1., &status);
    for (block = 0; block < NBLOCKS; block++)
    {
        v = VALUE(rank, step, block);
        log("  Write block %d, value %d to %s\n", block, v, FILENAME);
        set_vars(step, block);
        size_t start[2] = {offs1, offs2};
        for (i = 0; i < NVARS; i++)
        {
            adios2_set_selection(varW[i], 2, start, count);
            adios2_put(engineW, varW[i], a2, adios2_mode_sync);
        }
    }
    adios2_end_step(engineW);

    te = MPI_Wtime();
    if (rank == 0)
    {
        log("  Write time for step %d was %6.3lf seconds\n", step, te - tb);
    }
    MPI_Barrier(comm);
    return 0;
}

#define CHECK_VARINFO(VARNAME, NDIM, NSTEPS)                                                       \
    vi = adios2_inquire_variable(ioR, VARNAME);                                                    \
    if (vi == NULL)                                                                                \
    {                                                                                              \
        printE("No such variable: %s\n", VARNAME);                                                 \
        err = 101;                                                                                 \
        goto endread;                                                                              \
    }                                                                                              \
    size_t ndims;                                                                                  \
    adios2_variable_ndims(&ndims, vi);                                                             \
                                                                                                   \
    if (ndims != NDIM)                                                                             \
    {                                                                                              \
        printE("Variable %s has %zu dimensions, but expected %u\n", VARNAME, ndims, NDIM);         \
        err = 102;                                                                                 \
        goto endread;                                                                              \
    }                                                                                              \
    size_t steps;                                                                                  \
    adios2_variable_steps(&steps, vi);                                                             \
    if (steps != NSTEPS)                                                                           \
    {                                                                                              \
        printE("Variable %s has %zu steps, but expected %u\n", VARNAME, steps, NSTEPS);            \
        err = 103;                                                                                 \
        /*goto endread; */                                                                         \
    }

#define CHECK_SCALAR(VARNAME, VAR, VALUE, STEP)                                                    \
    if (VAR != VALUE)                                                                              \
    {                                                                                              \
        printE(#VARNAME " step %d: wrote %d but read %d\n", STEP, VALUE, VAR);                     \
        err = 104;                                                                                 \
        /*goto endread;*/                                                                          \
    }

#define CHECK_ARRAY(VARNAME, A, N, VALUE, STEP, BLOCK, i)                                          \
    for (i = 0; i < N; i++)                                                                        \
        if (A[i] != VALUE)                                                                         \
        {                                                                                          \
            printE("%s[%d] step %d block %d: wrote %d but read %d\n", VARNAME, i, STEP, BLOCK,     \
                   VALUE, A[i]);                                                                   \
            err = 104;                                                                             \
            /*goto endread;*/                                                                      \
            break;                                                                                 \
        }

void reset_readvars()
{
    int n;

    n = ldim1 * ldim2;
    memset(r2, -1, n * sizeof(int));
}

int read_file()
{
    adios2_variable *vi;
    int err = 0, v;
    int block, step, i; // loop variables
    int iMacro;         // loop variable in macros
    double tb, te;
    double tsb, ts; // time for just scheduling for one step/block

    size_t start[2] = {offs1, offs2};
    size_t count[2] = {ldim1, ldim2};

    reset_readvars();

    log("Read and check data in %s\n", FILENAME);
    adios2_engine *engineR = adios2_open(ioR, FILENAME, adios2_mode_read);
    if (engineR == NULL)
    {
        printE("Error at opening file %s for reading\n", FILENAME);
        return 1;
    }

    adios2_step_status status;
    adios2_begin_step(engineR, adios2_step_mode_read, -1., &status);

    log("  Check variable definitions... %s\n", FILENAME);
    tb = MPI_Wtime();
    for (i = 0; i < NVARS; i++)
    {
        CHECK_VARINFO(varnames[i], 2, NSTEPS)
    }
    MPI_Barrier(comm);
    te = MPI_Wtime();
    if (rank == 0)
    {
        log("  Time to check all vars' info: %6.3lf seconds\n", te - tb);
    }

    log("  Check variable content...\n");
    for (step = 0; step < NSTEPS; step++)
    {
        tb = MPI_Wtime();
        ts = 0;
        if (step > 0)
        {
            adios2_begin_step(engineR, adios2_step_mode_read, -1., &status);
        }
        for (block = 0; block < NBLOCKS; block++)
        {

            if (status == adios2_step_status_ok)
            {
                v = VALUE(rank, step, block);
                set_offsets(block);
                start[0] = offs1;
                start[1] = offs2;

                log("    Step %d block %d: value=%d\n", step, block, v);

                for (i = 0; i < NVARS; i++)
                {
                    tsb = MPI_Wtime();
                    adios2_variable *varH = adios2_inquire_variable(ioR, varnames[i]);
                    adios2_set_selection(varH, 2, start, count);
                    adios2_get(engineR, varH, r2, adios2_mode_sync);
                    ts += MPI_Wtime() - tsb;
                    CHECK_ARRAY(varnames[i], r2, ldim1 * ldim2, v, step, block, iMacro)
                }
            }
            else
            {
                printf("-- ERROR: Could not get Step %d, status = %d\n", i, status);
            }
        }
        adios2_end_step(engineR);
        MPI_Barrier(comm);
        te = MPI_Wtime();
        if (rank == 0)
        {
            log("  Read time for step %d was %6.3lfs\n", step, ts);
        }
    }

endread:

    adios2_close(engineR);
    MPI_Barrier(comm);
    return err;
}
