/*
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/*=================================================================
 * adiosreadc.c - read a variable/attribute from an ADIOS file
 *
 * Input: File, Group, Path, Starts, Counts, StepStart, StepCount, Verbose
 *    File:      uint64 (handler)
 *    Group:     uint64 (handler)
 *    Path:      string (variable name with full path)
 *    Starts:    int64 array
 *    Counts:    int64 array
 *    StepStart: int64
 *    StepCount: int64
 *    Verbose:   numeric (double)
 * Output: The variable/attribute as mxArray
 *
 *
 * Date: 2018/09/07
 * Author: Norbert Podhorszki <pnorbert@ornl.gov>
 *=================================================================*/

#include "adios2_c.h"
#include "matrix.h"
#include "mex.h"
#include <stdint.h> /* uint64_t and other int types */
#include <string.h> /* strlen */

static int verbose = 0;

mxArray *readdata(adios2_engine *fp, adios2_io *group, const char *path, mwSize in_noffsets,
                  const int64_t *in_offsets, const int64_t *in_counts, const int64_t in_stepstart,
                  const int64_t in_stepcount);
void errorCheck(int nlhs, int nrhs, const mxArray *prhs[]);
void checkDimSize(const int ndims, const size_t *dims);
char *getString(const mxArray *mxstr);
mxClassID adiostypeToMatlabClass(int adiostype, mxComplexity *complexity);
mxArray *createMatlabArray(int adiostype, size_t ndim, size_t *dims);
void recalc_offsets(const size_t ndim, const size_t *dims, mwSize in_noffsets,
                    const int64_t *in_offsets, const int64_t *in_counts, size_t *offsets,
                    size_t *counts);
void recalc_steps(const size_t varStepStart, const size_t varStepCount, const int64_t in_stepstart,
                  const int64_t in_stepcount, size_t *start, size_t *count);
static void swap_order(size_t n, size_t *array);
void printArrayInt64(size_t nelems, void *array);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char *fname;         /* file name */
    char *path;          /* name of variable */
    mwSize in_noffsets;  /* number of offsets provided, = number of counts */
    int64_t *in_offsets; /* offset array provided */
    int64_t *in_counts;  /* count array provided */
    char msg[512];       /* error messages from function calls */

    uint64_t fh, gh, *uint64p; /* file and group handlers and temp pointers */
    int32_t *int32p;

    adios2_io *group;  /* IO group object pointer */
    adios2_engine *fp; /* File handler pointer  */

    mxArray *out; /* output data array */

    errorCheck(nlhs, nrhs, prhs);

    /*mexPrintf("nrhs=%d  nlhs=%d\n", nrhs, nlhs);*/

    /***********************/
    /* 0. get verbose parameter first */
    verbose = (int)mxGetScalar(prhs[7]);
    if (verbose)
        mexPrintf("Verbosity level: %d\n", verbose);

    /* 1. get file handler */
    uint64p = (uint64_t *)mxGetData(prhs[0]);
    fp = (adios2_engine *)*uint64p;
    if (verbose)
        mexPrintf("File handler: \"%llu\"\n", (uint64_t)fp);

    /* 2. get IO group object handler */
    uint64p = (uint64_t *)mxGetData(prhs[1]);
    group = (adios2_io *)*uint64p;
    if (verbose)
        mexPrintf("IO group handler: \"%llu\"\n", (uint64_t)group);

    /*****************************************************************************************/
    /* 3. get other arguments: char variable name, int64 in_offsets[], int64
     * in_counts[] */
    path = getString(prhs[2]);
    if (verbose)
        mexPrintf("Variable name: \"%s\"\n", path);
    in_noffsets = mxGetNumberOfElements(prhs[3]);
    in_offsets = (int64_t *)mxGetData(prhs[3]);
    in_counts = (int64_t *)mxGetData(prhs[4]);

    int64_t in_stepstart = *(int64_t *)mxGetData(prhs[5]);
    int64_t in_stepcount = *(int64_t *)mxGetData(prhs[6]);
    if (verbose)
        mexPrintf("Step start: %lld\n", in_stepstart);
    if (verbose)
        mexPrintf("Step count: %lld\n", in_stepcount);

    /*********************************************************************/
    /* 5. read in variable */
    out = readdata(fp, group, path, in_noffsets, in_offsets, in_counts, in_stepstart, in_stepcount);
    if (nlhs >= 1)
    {
        plhs[0] = out;
    }

    mxFree(path);
}

mxArray *readdata(adios2_engine *fp, adios2_io *group, const char *path, mwSize in_noffsets,
                  const int64_t *in_offsets, const int64_t *in_counts, const int64_t in_stepstart,
                  const int64_t in_stepcount)
{
    /* FIXME: does not work for COMPLEX data because
       adios stores real and imaginary together as Fortran Complex type does,
       while matlab stores them as two separate float/double arrays
    */
    adios2_variable *avar; /* Variable object */
    size_t read_bytes;
    size_t arraysize;
    int i;
    char msg[512];
    adios2_variable *vinfo;

    /* get variable object */
    avar = adios2_inquire_variable(group, path);
    if (!avar)
    {
        mexErrMsgIdAndTxt("MATLAB:adiosreadc:variable %s does not exist", path);
    }

    size_t varStepStart;
    adios2_variable_steps_start(&varStepStart, avar);

    size_t varStepCount;
    adios2_variable_steps(&varStepCount, avar);

    size_t varNdim;
    adios2_variable_ndims(&varNdim, avar);

    size_t *varDim = (size_t *)malloc(varNdim * sizeof(size_t));
    adios2_variable_shape(varDim, avar);

    adios2_type adiostype;
    adios2_variable_type(&adiostype, avar);

    /* if (verbose) {
        mexPrintf("Var=%s type=%s C dimensions %dD [", path,
    adios_type_to_string(vinfo->type), vinfo->varNdim);
        for (i=0; i<vinfo->varNdim; i++)
            mexPrintf("%lld ", vinfo->varDim[i]);
        mexPrintf("]\n");
    }*/

    size_t mxndim = varNdim;
    if (varStepCount > 1)
    {
        ++mxndim;
    }
    size_t *mxdims = (size_t *)mxCalloc(mxndim, sizeof(size_t));
    for (i = 0; i < varNdim; i++)
    {
        mxdims[i] = varDim[i];
    }

    /* Flip dimensions from ADIOS-read-api/C/row-major order to
     * Matlab/Fortran/column-major order */
    swap_order(varNdim,
               mxdims); // swap only real varDim, not the steps if present

    /* extend offsets/counts if needed, change from 1..N to 0..N-1 indexing and
       interpret negative indices */
    size_t qoffsets[16], qcounts[16]; /* extended offsets/counts */
    size_t qstepstart = varStepStart; // default value in case there are no
                                      // steps for the variable
    size_t qstepcount = varStepCount;
    recalc_offsets(varNdim, mxdims, in_noffsets, in_offsets, in_counts, qoffsets, qcounts);
    recalc_steps(varStepStart, varStepCount, in_stepstart, in_stepcount, &qstepstart, &qstepcount);

    if (mxndim > varNdim)
    {
        if (verbose)
            mexPrintf("Add steps as extra dimension, start=%zu  count=%zu ", qstepstart,
                      qstepcount);
        qoffsets[varNdim] = qstepstart;
        qcounts[varNdim] = qstepcount; // steps become slowest dimension
    }

    checkDimSize(varNdim, qcounts);

    /* create Matlab array with the appropriate type and size */
    void *data;
    mxArray *out;
    if (adiostype == adios2_type_string)
    {
        /* Matlab string != char array of C, so handle separately */
        data = (void *)mxCalloc(65536, sizeof(char));
        mexPrintf("Create C string reading in a string with ptr %p\n", data);
    }
    else
    {
        if (verbose)
        {
            mexPrintf("Create %d-D Matlab array ", mxndim);
            printArrayInt64(mxndim, qcounts);
            mexPrintf("\n");
        }
        out = createMatlabArray(adiostype, mxndim, qcounts);
        data = (void *)mxGetData(out);
    }

    /* Flip offsets/counts from Matlab/Fortran/column-major order to
     * ADIOS-read-api/C/row-major order */
    swap_order(varNdim, qoffsets); // again, leave out the steps from the swap
    swap_order(varNdim, qcounts);

    if (verbose)
    {
        mexPrintf("Set selection for variable: start = ");
        printArrayInt64(varNdim, qoffsets);
        mexPrintf("  count = ");
        printArrayInt64(varNdim, qcounts);
        mexPrintf("\n");
        mexPrintf("Set step-selection for variable: start = %zu  count = %zu\n", qstepstart,
                  qstepcount);
    }

    if (varNdim > 0)
    {
        adios2_set_selection(avar, varNdim, qoffsets, qcounts);
    }
    adios2_set_step_selection(avar, qstepstart, qstepcount);

    /* read in data */
    if (verbose)
        mexPrintf("Read in data\n");

    adios2_get(fp, avar, data, adios2_mode_sync);

    if (adiostype == adios2_type_string)
    {
        mexPrintf("Create Matlab string from C string [%s]\n", (char *)data);
        out = mxCreateString((char *)data);
        mxFree(data);
    }

    free(varDim);

    return out;
}

void errorCheck(int nlhs, int nrhs, const mxArray *prhs[])
{
    /* Assume that we are called from adiosread.m which checks the arguments
     * already */
    /* Check for proper number of arguments. */

    if (nrhs != 8)
    {
        mexErrMsgIdAndTxt("MATLAB:adiosreadc:rhs", "This function needs exactly 6 arguments: File, "
                                                   "Group, Varpath, Offsets, Counts, StepStart, "
                                                   "StepCount, Verbose");
    }

    if (!mxIsUint64(prhs[0]))
    {
        mexErrMsgIdAndTxt("MATLAB:adiosreadc:rhs", "First arg must be an uint64 handler.");
    }

    if (!mxIsUint64(prhs[1]))
    {
        mexErrMsgIdAndTxt("MATLAB:adiosreadc:rhs", "Second arg must be an uint64 handler.");
    }

    if (!mxIsChar(prhs[2]))
    {
        mexErrMsgIdAndTxt("MATLAB:adiosreadc:rhs", "Third arg must be a string.");
    }

    if (nlhs > 1)
    {
        mexErrMsgIdAndTxt("MATLAB:adiosreadc:lhs", "Too many output arguments.");
    }
}

void checkDimSize(const int ndims, const size_t *dims)
{
    /* Make sure that it is safe to cast dim to mwSize (MWSIZE_MAX is limited)*/
    for (int i = 0; i < ndims; ++i)
    {
        if (dims[i] > MWSIZE_MAX)
        {
            mexErrMsgIdAndTxt("MATLAB:adiosreadc:dimensionTooLarge",
                              "The selected dimension size, %zu, is larger than the "
                              "maximum supported value of mwSize, %u, so we cannot create "
                              "the result array\n",
                              dims[i], MWSIZE_MAX);
        }
    }
}

/** Make a C char* string from a Matlab string */
char *getString(const mxArray *mxstr)
{
    mwSize buflen;
    char *str;
    /* Allocate enough memory to hold the converted string. */
    buflen = mxGetNumberOfElements(mxstr) + 1;
    str = mxCalloc(buflen, sizeof(char));
    /* Copy the string data from string_array_ptr and place it into buf. */
    if (mxGetString(mxstr, str, buflen) != 0)
        mexErrMsgTxt("Could not convert string data from the file name.");
    return str;
}

/** return the appropriate class for an adios type (and complexity too) */
mxClassID adiostypeToMatlabClass(adios2_type adiostype, mxComplexity *complexity)
{
    *complexity = mxREAL;
    switch (adiostype)
    {
    case adios2_type_uint8_t:
        return mxUINT8_CLASS;
    case adios2_type_int8_t:
        return mxINT8_CLASS;

    case adios2_type_string:
        /* case adios2_type_string_array: */
        return mxCHAR_CLASS;

    case adios2_type_uint16_t:
        return mxUINT16_CLASS;
    case adios2_type_int16_t:
        return mxINT16_CLASS;

    case adios2_type_uint32_t:
        return mxUINT32_CLASS;
    case adios2_type_int32_t:
        return mxINT32_CLASS;

    case adios2_type_uint64_t:
        return mxUINT64_CLASS;
    case adios2_type_int64_t:
        return mxINT64_CLASS;

    case adios2_type_float:
        return mxSINGLE_CLASS;
    case adios2_type_double:
        return mxDOUBLE_CLASS;

    case adios2_type_float_complex: /* 8 bytes */
        *complexity = mxCOMPLEX;
        return mxSINGLE_CLASS;
    case adios2_type_double_complex: /*  16 bytes */
        *complexity = mxCOMPLEX;
        return mxDOUBLE_CLASS;

    default:
        mexErrMsgIdAndTxt("MATLAB:adiosopenc.c:dimensionTooLarge",
                          "Adios type id=%d not supported in matlab.\n", adiostype);
        break;
    }
    return 0; /* just to avoid warnings. never executed */
}

/* create an N-dim array with given type and dimensions */
mxArray *createMatlabArray(int adiostype, size_t ndim, size_t *dims)
{
    mxClassID mxtype;
    mxArray *arr;
    mwSize mNdim;
    mwSize mDims[16];
    mxComplexity ComplexFlag;
    int i;
    int isTypeComplex;

    /* convert ints to mwSizes */
    if (ndim > 0)
    {
        mNdim = (mwSize)ndim;
        for (i = 0; i < ndim; i++)
            mDims[i] = (mwSize)dims[i];
    }
    else
    {
        /* 0 dim: scalar value -> 1-by-1 Matlab array */
        mNdim = 2;
        mDims[0] = 1;
        mDims[1] = 1;
    }

    /* get type */
    mxtype = adiostypeToMatlabClass(adiostype, &ComplexFlag);

    /* create array */
    arr = mxCreateNumericArray(mNdim, mDims, mxtype, ComplexFlag);

    if (verbose)
        mexPrintf("Array for adios type %d is created\n", adiostype);

    return arr;
}

/** - extend offset/count arrays to ndims if needed
    - recalculate "from 1" Matlab indices to "from 0" C indices
    - recalculate the negative indices
    !!! Provide the output arrays in the caller !!!
*/
void recalc_offsets(const size_t ndim, const size_t *dims, mwSize in_noffsets,
                    const int64_t *in_offsets, const int64_t *in_counts, size_t *offsets,
                    size_t *counts)
{
    int i;
    for (i = 0; i < ndim; i++)
    {
        if ((mwSize)i < in_noffsets)
        {
            if (in_offsets[i] < 0) /* negative offset means last-|offset| */
                offsets[i] = dims[i] + in_offsets[i];
            else
                offsets[i] = in_offsets[i] - 1; /* C index start from 0 */

            if (in_counts[i] < 0) /* negative count means last-|count|+1-start */
                counts[i] = dims[i] + in_counts[i] - offsets[i] + 1;
            else
                counts[i] = in_counts[i];
        }
        else
        {
            /* extend offset/count array to match variable's dimensions */
            if (verbose)
                mexPrintf("Extend offset/counts for dim %d: offset=%d count=%d\n", i, 0, dims[i]);
            offsets[i] = 0;
            counts[i] = dims[i];
        }
    }
}

void recalc_steps(const size_t varStepStart, const size_t varStepCount, const int64_t in_stepstart,
                  const int64_t in_stepcount, size_t *start, size_t *count)
{
    /* handle steps for variables with multiple steps */
    if (in_stepstart < 0) /* negative offset means last step -|start| */
        *start = varStepCount + in_stepstart;
    else
        *start = in_stepstart - 1; /* C index start from 0 */

    if (in_stepcount < 0) /* negative offset means last step -|start| */
    {
        *count = varStepCount + in_stepcount + 1 - *start;
    }
    else
        *count = in_stepcount; /* C index start from 0 */
}

/* Reverse the order in an array in place.
   use swapping from Matlab/Fortran/column-major order to
   ADIOS-read-api/C/row-major order and back
*/
static void swap_order(size_t n, size_t *array)
{
    int i, tmp;
    for (i = 0; i < n / 2; i++)
    {
        tmp = array[i];
        array[i] = array[n - 1 - i];
        array[n - 1 - i] = tmp;
    }
}

void printArrayInt64(size_t nelems, void *array)
{
    int64_t *a = (int64_t *)array;
    size_t i;
    mexPrintf("[");
    for (i = 0; i < nelems; i++)
        mexPrintf("%lld%s", a[i], (i < nelems - 1 ? " " : ""));
    mexPrintf("]");
}
