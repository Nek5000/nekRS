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

#include <chrono>
#include <vector>

#include <gtest/gtest.h>

#include <adios2/common/ADIOSConfig.h>

#if ADIOS2_USE_MPI
#include <mpi.h>
#endif

#include "adios2_c.h"
#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef DMALLOC
#include "dmalloc.h"
#endif

struct RunParams
{
    size_t nvars;
    size_t nblocks;
    size_t nsteps;
    RunParams(size_t nv, size_t nb, size_t ns) : nvars{nv}, nblocks{nb}, nsteps{ns} {};
};

/* This function is executed by INSTANTIATE_TEST_SUITE_P
   before main() and MPI_Init()!!! */
std::vector<RunParams> CreateRunParams()
{
    std::vector<RunParams> params;
    // 1 variable
    params.push_back(RunParams(1, 1, 1));
    params.push_back(RunParams(1, 1, 2));
    params.push_back(RunParams(1, 2, 1));
    params.push_back(RunParams(1, 2, 2));
    // 2 variables
    params.push_back(RunParams(2, 1, 1));
    params.push_back(RunParams(2, 1, 2));
    params.push_back(RunParams(2, 2, 1));
    params.push_back(RunParams(2, 2, 2));

    // These pass
    params.push_back(RunParams(5, 14, 1));
    params.push_back(RunParams(8, 8, 1));

    // These fail
    // params.push_back(RunParams(5, 15, 1));
    // params.push_back(RunParams(8, 9, 1));

    return params;
}

#define log(...)                                                                                   \
    fprintf(stderr, "[rank=%3.3d, line %d]: ", rank, __LINE__);                                    \
    fprintf(stderr, __VA_ARGS__);                                                                  \
    fflush(stderr);
#define printE(...)                                                                                \
    fprintf(stderr, "[rank=%3.3d, line %d]: ERROR: ", rank, __LINE__);                             \
    fprintf(stderr, __VA_ARGS__);                                                                  \
    fflush(stderr);

#define VALUE(rank, step, block) (step * 10000 + 10 * rank + block)

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
    if (ndims != NDIM)                                                                             \
    {                                                                                              \
        printE("Variable %s has %zu dimensions, but expected %u\n", VARNAME, ndims, NDIM);         \
        err = 102;                                                                                 \
        goto endread;                                                                              \
    }

#define CHECK_SCALAR(VARNAME, VAR, VALUE, STEP)                                                    \
    if (VAR != VALUE)                                                                              \
    {                                                                                              \
        printE(#VARNAME " step %zu: wrote %zu but read %zu\n", STEP, VALUE, VAR);                  \
        err = 104;                                                                                 \
        /*goto endread;*/                                                                          \
    }

#define CHECK_ARRAY(VARNAME, A, N, VALUE, STEP, BLOCK)                                             \
    for (int i = 0; i < static_cast<int>(N); i++)                                                  \
        if (A[i] != VALUE)                                                                         \
        {                                                                                          \
            printE("%s[%d] step %d block %zu: wrote %d but read %d\n", VARNAME, i, STEP, BLOCK,    \
                   VALUE, A[i]);                                                                   \
            err = 104;                                                                             \
            /*goto endread;*/                                                                      \
            break;                                                                                 \
        }

#if ADIOS2_USE_MPI
MPI_Comm comm = MPI_COMM_WORLD;
#endif
int rank = 0;
int numprocs = 1;

class TestManyVars : public ::testing::TestWithParam<RunParams>
{
public:
    TestManyVars() = default;
    ~TestManyVars() = default;

    size_t NVARS = 1;
    size_t NBLOCKS = 1;
    size_t NSTEPS = 1;
    int REDEFINE = 0; // 1: delete and redefine variable definitions at each
                      // step to test adios_delete_vardefs()
    char FILENAME[256];

    /* Variables to write */
    int *a2 = nullptr;

    static const size_t ldim1 = 5;
    static const size_t ldim2 = 5;
    size_t gdim1 = 0, gdim2 = 0;
    size_t offs1 = 0, offs2 = 0;

    /* ADIOS2 variables to be accessed from multiple functions */
    adios2_io *ioW, *ioR;
    adios2_engine *engineW;
    adios2_variable **varW; // array to hold all variable definitions

    /* Variables to read */
    int *r2;
    char **varnames;

    void alloc_vars()
    {
        size_t n = ldim1 * ldim2;
        a2 = (int *)malloc(n * sizeof(int));
        r2 = (int *)malloc(n * sizeof(int));
        varW = (adios2_variable **)malloc(NVARS * sizeof(adios2_variable *));

        varnames = (char **)malloc(NVARS * sizeof(char *));
        for (size_t i = 0; i < NVARS; i++)
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

        char fmt[32];
        snprintf(fmt, sizeof(fmt), "v%%%d.%dd", digit, digit);
        for (size_t i = 0; i < NVARS; i++)
        {
            snprintf(varnames[i], 16, fmt, i);
        }
    }

    void set_gdim()
    {
        gdim1 = numprocs * ldim1;
        gdim2 = NBLOCKS * ldim2;
    }

    void set_offsets(size_t block)
    {
        offs1 = rank * ldim1;
        offs2 = block * ldim2;
    }

    void set_vars(int step, int block)
    {
        int v = VALUE(rank, step, block);

        set_offsets(block);

        size_t n = ldim1 * ldim2;
        for (size_t i = 0; i < n; i++)
            a2[i] = v;
    }

    void fini_vars()
    {
        free(a2);
        free(r2);
        free(varW);
        for (size_t i = 0; i < NVARS; i++)
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

    int Test(RunParams p, bool redefineVars)
    {
        int err;

        NVARS = p.nvars;
        NBLOCKS = p.nblocks;
        NSTEPS = p.nsteps;
        REDEFINE = redefineVars;
#if ADIOS2_USE_MPI
        snprintf(FILENAME, sizeof(FILENAME), "manyVars.%zu_%zu_%zu%s_MPI.bp", NVARS, NBLOCKS,
                 NSTEPS, REDEFINE ? "_redefine" : "");
#else
        snprintf(FILENAME, sizeof(FILENAME), "manyVars.%zu_%zu_%zu%s.bp", NVARS, NBLOCKS, NSTEPS,
                 REDEFINE ? "_redefine" : "");
#endif

        alloc_vars();
#if ADIOS2_USE_MPI
        adios2_adios *adiosH = adios2_init_mpi(MPI_COMM_WORLD);
#else
        adios2_adios *adiosH = adios2_init_serial();
#endif

        ioW = adios2_declare_io(adiosH, "multiblockwrite"); // group for writing
        ioR = adios2_declare_io(adiosH, "multiblockread");  // group for reading
        set_gdim();

        if (rank == 0)
        {
            log("Test %zu Variables, %zu Blocks, %zu Steps\n", NVARS, NBLOCKS, NSTEPS);
        }

        engineW = adios2_open(ioW, FILENAME, adios2_mode_write);

        err = 0;
        for (size_t i = 0; i < NSTEPS; i++)
        {
            if (!err)
            {
                if (i == 0 || REDEFINE)
                {
                    printf("-- Define variables.\n");
                    define_vars();
                }

                err = write_file(i);

                std::cout << "After writem err is " << err << std::endl;
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

        std::cout << "After read err is " << err << std::endl;
        adios2_finalize(adiosH);
        fini_vars();
        return err;
    }

    void define_vars()
    {
        size_t shape[2] = {gdim1, gdim2};
        size_t count[2] = {ldim1, ldim2};
        size_t start[2] = {0, 0};

        /* One variable definition for many blocks.
         * Offsets will change at writing for each block. */
        for (size_t i = 0; i < NVARS; i++)
        {
            varW[i] = adios2_define_variable(ioW, varnames[i], adios2_type_int32_t, 2, shape, start,
                                             count, adios2_constant_dims_false);
        }
    }

    int write_file(int step)
    {
        int v;
        size_t count[2] = {ldim1, ldim2};

        log("Write step %d to %s\n", step, FILENAME);
        auto tb = std::chrono::high_resolution_clock::now();

        adios2_step_status status;
        adios2_begin_step(engineW, adios2_step_mode_append, -1., &status);
        for (size_t block = 0; block < NBLOCKS; block++)
        {
            v = VALUE(rank, step, block);
            log("  Write block %d, value %d to %s\n", static_cast<int>(block), v, FILENAME);
            set_vars(step, block);
            size_t start[2] = {offs1, offs2};
            for (size_t i = 0; i < NVARS; i++)
            {
                adios2_set_selection(varW[i], 2, start, count);
                adios2_put(engineW, varW[i], a2, adios2_mode_sync);
            }
        }
        adios2_perform_puts(engineW);
        adios2_end_step(engineW);

        auto te = std::chrono::high_resolution_clock::now();
        if (rank == 0)
        {
            log("  Write time for step %d was %6.3lf seconds\n", step,
                double(std::chrono::duration_cast<std::chrono::seconds>(te - tb).count()));
        }
#if ADIOS2_USE_MPI
        MPI_Barrier(comm);
#endif
        return 0;
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
        int err = 0;
        std::chrono::time_point<std::chrono::high_resolution_clock> tb, te;
        std::chrono::duration<double> ts; // time for just scheduling for one step/block

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
        tb = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < NVARS; i++)
        {
            CHECK_VARINFO(varnames[i], 2, NSTEPS)
        }
#if ADIOS2_USE_MPI
        MPI_Barrier(comm);
#endif
        te = std::chrono::high_resolution_clock::now();
        if (rank == 0)
        {
            log("  Time to check all vars' info: %6.3lf seconds\n",
                double(std::chrono::duration_cast<std::chrono::seconds>(te - tb).count()));
        }

        log("  Check variable content...\n");
        for (int step = 0; step < static_cast<int>(NSTEPS); step++)
        {
            tb = std::chrono::high_resolution_clock::now();
            ts = std::chrono::duration<double>::zero();
            if (step > 0)
            {
                adios2_begin_step(engineR, adios2_step_mode_read, -1., &status);
            }
            for (size_t block = 0; block < NBLOCKS; block++)
            {
                if (status == adios2_step_status_ok)
                {
                    int v = VALUE(rank, step, block);
                    set_offsets(block);
                    start[0] = offs1;
                    start[1] = offs2;

                    log("    Step %d block %zu: value=%d\n", step, block, v);

                    for (size_t i = 0; i < NVARS; i++)
                    {
                        auto tsb = std::chrono::high_resolution_clock::now();
                        adios2_variable *varH = adios2_inquire_variable(ioR, varnames[i]);
                        adios2_set_selection(varH, 2, start, count);
                        adios2_get(engineR, varH, r2, adios2_mode_sync);
                        auto tse = std::chrono::high_resolution_clock::now();
                        ts += tse - tsb;
                        CHECK_ARRAY(varnames[i], r2, ldim1 * ldim2, v, step, block)
                    }
                }
                else
                {
                    printf("-- ERROR: Could not get Step %d, status = %d\n", step, status);
                }
            }
            adios2_end_step(engineR);
#if ADIOS2_USE_MPI
            MPI_Barrier(comm);
#endif
            te = std::chrono::high_resolution_clock::now();
            if (rank == 0)
            {
                log("  Read time for step %d was %6.3lfs\n", step,
                    double(std::chrono::duration_cast<std::chrono::seconds>(ts).count()));
            }
        }

    endread:

        adios2_close(engineR);
#if ADIOS2_USE_MPI
        MPI_Barrier(comm);
#endif
        return err;
    }
};

TEST_P(TestManyVars, DontRedefineVars)
{
    RunParams p = GetParam();
    int err = Test(p, false);
    ASSERT_EQ(err, 0);
}

INSTANTIATE_TEST_SUITE_P(NxM, TestManyVars, ::testing::ValuesIn(CreateRunParams()));

//******************************************************************************
// main
//******************************************************************************

int main(int argc, char **argv)
{
#if ADIOS2_USE_MPI
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
#endif
    ::testing::InitGoogleTest(&argc, argv);

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
#endif

    int result;
    result = RUN_ALL_TESTS();

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif
    return result;
}
