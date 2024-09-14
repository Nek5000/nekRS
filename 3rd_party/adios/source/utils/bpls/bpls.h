/*
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include "adios2/common/ADIOSConfig.h"
#include "adios2/common/ADIOSMacros.h"
#include "adios2/common/ADIOSTypes.h"
#include "adios2/core/ADIOS.h"
#include "adios2/core/Engine.h"
#include "adios2/core/IO.h"
#include "adios2/core/Variable.h"
#include "adios2/helper/adiosFunctions.h"

#include "adios2/core/Info.h"

#include <map>
#include <string>

namespace adios2
{
namespace utils
{

/* definitions for bpls.c */
#define myfree(p)                                                                                  \
    if (p)                                                                                         \
    {                                                                                              \
        free(p);                                                                                   \
        p = NULL;                                                                                  \
    }

#define CUT_TO_BYTE(x) (x < 0 ? 0 : (x > 255 ? 255 : x))

#define MAX_DIMS 16
#define MAX_MASKS 10
#define MAX_BUFFERSIZE (10 * 1024 * 1024)

struct Entry
{
    DataType typeName;
    bool isVar;
    union
    {
        core::VariableBase *var;
        core::AttributeBase *attr;
    };
    Entry(DataType type, core::VariableBase *v) : typeName(type), isVar(true), var(v) {}
    Entry(DataType type, core::AttributeBase *a) : typeName(type), isVar(false), attr(a) {}
};

// how to print one data item of an array
// enum PrintDataType {STRING, INT, FLOAT, DOUBLE, COMPLEX};

char *mystrndup(const char *s, size_t n);
void init_globals();
void processDimSpecs();
void parseDimSpec(const std::string &str, int64_t *dims);
int parseAccuracy();
int compile_regexp_masks(void);
void printSettings(void);
int doList(std::string path);
void mergeLists(int nV, char **listV, int nA, char **listA, char **mlist, bool *isVar);

template <class T>
int printVariableInfo(core::Engine *fp, core::IO *io, core::Variable<T> *variable);

template <class T>
int readVar(core::Engine *fp, core::IO *io, core::Variable<T> *variable);

template <class T>
int readVarBlock(core::Engine *fp, core::IO *io, core::Variable<T> *variable, int blockid);

template <class T>
size_t relative_to_absolute_step(core::Engine *fp, core::Variable<T> *variable,
                                 const size_t relstep);
template <class T>
Dims get_global_array_signature(core::Engine *fp, core::IO *io, core::Variable<T> *variable);
template <class T>
std::pair<size_t, Dims> get_local_array_signature(core::Engine *fp, core::IO *io,
                                                  core::Variable<T> *variable);

int cmpstringp(const void *p1, const void *p2);
bool grpMatchesMask(char *name);
bool matchesAMask(const char *name);
int print_start(const std::string &fnamestr);
void print_slice_info(core::VariableBase *variable, bool timed, uint64_t *s, uint64_t *c,
                      Dims count);
int print_data(const void *data, int item, DataType adiosvartypes, bool allowformat,
               bool char_star_string = false);

/* s is a character array not necessarily null terminated.
 * return false on OK print, true if it not XML (not printed)*/
bool print_data_xml(const char *s, const size_t length);

int print_dataset(const void *data, const DataType vartype, uint64_t *s, uint64_t *c, int tdims,
                  int *ndigits);
void print_endline(void);
void print_stop(void);
int print_data_hist(core::VariableBase *vi, char *varname);
int print_data_characteristics(void *min, void *max, double *avg, double *std_dev,
                               DataType adiosvartypes, bool allowformat);

template <class T>
void print_decomp(core::Engine *fp, core::IO *io, core::Variable<T> *variable);

template <class T>
void print_decomp_singlestep(core::Engine *fp, core::IO *io, core::Variable<T> *variable);
// close namespace
}
}
