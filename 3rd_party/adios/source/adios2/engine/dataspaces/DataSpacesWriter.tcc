/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * DataSpacesWriter.tcc
 *
 *  Created on: Dec 5, 2018
 *      Author: Pradeep Subedi
 *				pradeep.subedi@rutgers.edu
 */
#ifndef ADIOS2_ENGINE_DATASPACES_DATASPACESWRITER_TCC_
#define ADIOS2_ENGINE_DATASPACES_DATASPACESWRITER_TCC_

#include <memory>

#include "DataSpacesWriter.h"
#include "adios2/helper/adiosFunctions.h"
#include "adios2/toolkit/dataspaces/DSpacesConfig.h"
#include "adios2/toolkit/dataspaces/ds.h"
#include "adios2/toolkit/dataspaces/ds_data.h"
#ifdef HAVE_DSPACES2
#include "dspaces.h"
#else
#include "dataspaces.h"
#endif /* HAVE_DSPACES2 */

namespace adios2
{
namespace core
{
namespace engine
{

template <class T>
void DataSpacesWriter::DoPutSyncCommon(Variable<T> &variable, const T *values)
{

    uint64_t lb_in[MAX_DS_NDIM], ub_in[MAX_DS_NDIM], gdims_in[MAX_DS_NDIM];

    std::vector<uint64_t> dims_vec;

    unsigned int version;
    version = m_CurrentStep;
    int ndims = std::max(variable.m_Shape.size(), variable.m_Count.size());
    bool isOrderC = m_IO.m_ArrayOrder == adios2::ArrayOrdering::RowMajor;
    /* Order of dimensions: in DataSpaces: fast --> slow --> slowest
           For example:
           Fortran: i,j,k --> i, j, k  = lb[0], lb[1], lb[2]
                    i,j   --> i, j     = lb[0], lb[1]
                    i     --> i        = lb[0]
           C:       i,j,k --> k, j, i  = lb[2], lb[1], lb[0]
                    i,j   --> j, i     = lb[1], lb[0]
                    i     --> i        = lb[0]
        */
    if (variable.m_SingleValue)
    {
#ifdef HAVE_DSPACES2
        gdims_in[0] = 1;
#else
        gdims_in[0] = dspaces_get_num_space_server();
#endif /* HAVE_DSPACES2 */
        lb_in[0] = 0;
        ub_in[0] = 0;
        ndims = 1;
        dims_vec.push_back(0);
        ndim_vector.push_back(0);
    }
    else
    {
        ndim_vector.push_back(ndims);
        if (isOrderC)
        {
            for (int i = 0; i < ndims; i++)
            {
                gdims_in[i] = static_cast<uint64_t>(variable.m_Shape[ndims - i - 1]);
                dims_vec.push_back(gdims_in[i]);
                lb_in[i] = static_cast<uint64_t>(variable.m_Start[ndims - i - 1]);
                ub_in[i] = static_cast<uint64_t>(variable.m_Start[ndims - i - 1] +
                                                 variable.m_Count[ndims - i - 1] - 1);
            }
        }
        else
        {

            for (int i = 0; i < ndims; i++)
            {
                gdims_in[i] = static_cast<uint64_t>(variable.m_Shape[i]);
                dims_vec.push_back(gdims_in[i]);
                lb_in[i] = static_cast<uint64_t>(variable.m_Start[i]);
                ub_in[i] = static_cast<uint64_t>(variable.m_Start[i] + variable.m_Count[i] - 1);
            }
        }
    }
    gdims_vector.push_back(dims_vec);
    int varType;
    auto itType = varType_to_ds.find(ToString(variable.m_Type));
    if (itType == varType_to_ds.end())
    {
        varType = 2;
        fprintf(stderr, "variable Type not found. Using Integer as data type");
        // Might have to fix for complex data types
    }
    else
    {
        varType = itType->second;
    }
    elemSize_vector.push_back(varType);
    std::string ds_in_name = f_Name;
    ds_in_name += variable.m_Name;
    v_name_vector.push_back(variable.m_Name);
    char *var_str = new char[ds_in_name.length() + 1];
    strcpy(var_str, ds_in_name.c_str());
    variable.SetData(values);

    std::string l_Name = ds_in_name + std::to_string(version);
    char *cstr = new char[l_Name.length() + 1];
    strcpy(cstr, l_Name.c_str());

#ifdef HAVE_DSPACES2
    dspaces_client_t *client = get_client_handle();
    dspaces_define_gdim(*client, var_str, ndims, gdims_in);
    dspaces_put(*client, var_str, version, variable.m_ElementSize, ndims, lb_in, ub_in, values);
#else
    dspaces_define_gdim(var_str, ndims, gdims_in);
    dspaces_put(var_str, version, variable.m_ElementSize, ndims, lb_in, ub_in, values);
    dspaces_put_sync();
    dspaces_put_sync();
#endif /* HAVE_DSPACES2 */
    delete[] cstr;
    delete[] var_str;
}

} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif
