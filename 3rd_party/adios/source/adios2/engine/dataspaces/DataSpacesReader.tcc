/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * DataSpacesReader.tcc
 *
 *  Created on: Dec 5, 2018
 *      Author: Pradeep Subedi
 *				pradeep.subedi@rutgers.edu
 */
#ifndef ADIOS2_ENGINE_DATASPACES_DATASPACESREADER_TCC_
#define ADIOS2_ENGINE_DATASPACES_DATASPACESREADER_TCC_

#include <memory>

#include "DataSpacesReader.h"
#include "adios2/helper/adiosFunctions.h" //CSVToVector
#include "adios2/toolkit/dataspaces/DSpacesConfig.h"
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
void DataSpacesReader::AddVar(core::IO &io, std::string const &name, Dims shape)
{
    core::Variable<T> *v = io.InquireVariable<T>(name);
    if (NULL == v)
    {
        Dims zeros(shape.size(), 0);

        try
        {
            auto &foo = io.DefineVariable<T>(name, shape, zeros, shape);
        }
        catch (std::exception &e)
        {
            // invalid variable, do not define
            printf("WARNING: IO is not accepting definition of variable: %s. "
                   "Skipping. \n",
                   name.c_str());
        }
    }
    else
    {
        v->m_AvailableStepsCount++;
    }
}

template <class T>
void DataSpacesReader::ReadDsData(Variable<T> &variable, T *data, int version)
{
    uint64_t lb_in[MAX_DS_NDIM], ub_in[MAX_DS_NDIM], gdims_in[MAX_DS_NDIM];
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
    if (variable.m_Shape.size() == 0)
    {
#ifndef HAVE_DSPACES2
        gdims_in[0] = dspaces_get_num_space_server();
#else
        gdims_in[0] = 1;
#endif
        lb_in[0] = 0;
        ub_in[0] = 0;
        ndims = 1;
    }
    else
    {
        if (isOrderC)
        {
            for (int i = 0; i < ndims; i++)
            {
                gdims_in[i] = static_cast<uint64_t>(variable.m_Shape[ndims - i - 1]);
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
                lb_in[i] = static_cast<uint64_t>(variable.m_Start[i]);
                ub_in[i] = static_cast<uint64_t>(variable.m_Start[i] + variable.m_Count[i] - 1);
            }
        }
    }

    std::string ds_in_name = f_Name;
    ds_in_name += variable.m_Name;
    char *var_str = new char[ds_in_name.length() + 1];
    strcpy(var_str, ds_in_name.c_str());

    std::string l_Name = ds_in_name + std::to_string(version);
    char *cstr = new char[l_Name.length() + 1];
    strcpy(cstr, l_Name.c_str());
#ifdef HAVE_DSPACES2
    dspaces_client_t *client = get_client_handle();
    dspaces_define_gdim(*client, var_str, ndims, gdims_in);
    dspaces_get(*client, var_str, version, variable.m_ElementSize, ndims, lb_in, ub_in,
                (void *)data, -1);
#else
    dspaces_define_gdim(var_str, ndims, gdims_in);
    dspaces_get(var_str, version, variable.m_ElementSize, ndims, lb_in, ub_in, (void *)data);
#endif /* HAVE_DSPACES2 */
    delete[] cstr;
    delete[] var_str;
}

template <class T>
void DataSpacesReader::GetDeferredCommon(Variable<T> &variable, T *data)
{
    m_DeferredStack.push_back(variable.m_Name);
    variable.SetData(data);
}

template <class T>
void DataSpacesReader::GetSyncCommon(Variable<T> &variable, T *data)
{
    ReadDsData(variable, data, m_CurrentStep);
}

} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif
