/*
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/*=================================================================
 * adiosclosec.c - Close an ADIOS file
 *
 * Input: fh, ghs,  Verbose
 *    fh:       int64 adios file handler
 *    adiosobj: int64 adios object handler
 *    Verbose:  numeric (double)
 *
 * Date: 2018/09/07
 * Author: Norbert Podhorszki <pnorbert@ornl.gov>
 *=================================================================*/

#include "adios2_c.h"
#include "mex.h"
#include <stdint.h> /* uint64_t and other int types */
#include <string.h> /* memcpy */

static int verbose = 0;

mxClassID adiostypeToMatlabClass(int type, mxComplexity *complexity);
void errorCheck(int nlhs, int nrhs, const mxArray *prhs[]);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int status;
    char msg[512]; /* error messages from function calls */
    mwSize ngroups;
    uint64_t *uint64p; /* temp pointers */
    int i;
    adios2_adios *adiosobj; /* ADIOS object pointer */
    adios2_engine *fp;      /* File handler pointer  */

    errorCheck(nlhs, nrhs, prhs);

    /*mexPrintf("nrhs=%d  nlhs=%d\n", nrhs, nlhs);*/

    /***********************/
    /* 0. get verbose parameter first */
    verbose = (int)mxGetScalar(prhs[2]);
    if (verbose)
        mexPrintf("Verbosity level: %d\n", verbose);

    /* 1. get file handler */
    uint64p = (uint64_t *)mxGetData(prhs[0]);
    fp = (adios2_engine *)*uint64p;
    if (verbose)
        mexPrintf("File handler: \"%llu\"\n", *uint64p);

    /* 2. get adios object handler */
    uint64p = (uint64_t *)mxGetData(prhs[1]);
    adiosobj = (adios2_adios *)*uint64p;
    if (verbose)
        mexPrintf("ADIOS handler: \"%llu\"\n", *uint64p);

    size_t step;
    adios2_current_step(&step, fp);
    if (verbose)
        mexPrintf("Current step in file handler: %zu\n", step);

    if (verbose)
        mexPrintf("Close file handler\n");
    adios2_close(fp);

    if (verbose)
        mexPrintf("Destroy ADIOS object\n");
    adios2_finalize(adiosobj);

    if (verbose)
        mexPrintf("return from adiosclosec\n");
}

void errorCheck(int nlhs, int nrhs, const mxArray *prhs[])
{
    /* Assume that we are called from adiosread.m which checks the arguments
     * already */
    /* Check for proper number of arguments. */

    if (nrhs != 3)
    {
        mexErrMsgIdAndTxt("MATLAB:adiosclosec:rhs", "This function needs exactly 3 arguments: "
                                                    "FileHandler, ADIOSHandler, Verbose");
    }

    if (!mxIsUint64(prhs[0]))
    {
        mexErrMsgIdAndTxt("MATLAB:adiosclosec:rhs",
                          "First arg must be an uint64 handler to an ADIOS file .");
    }

    if (!mxIsUint64(prhs[1]))
    {
        mexErrMsgIdAndTxt("MATLAB:adiosclosec:rhs",
                          "Second arg must be an uint64 adios object handler.");
    }

    if (!mxIsNumeric(prhs[2]))
    {
        mexErrMsgIdAndTxt("MATLAB:adiosclosec:rhs", "Third arg must be a number.");
    }

    if (nlhs > 0)
    {
        mexErrMsgIdAndTxt("MATLAB:adiosclosec:lhs",
                          "This function does not have output arguments.");
    }
}
