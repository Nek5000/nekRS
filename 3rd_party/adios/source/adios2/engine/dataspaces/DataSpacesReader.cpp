/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * DataSpacesReader.cpp
 *
 *  Created on: Dec 5, 2018
 *      Author: Pradeep Subedi
 *      		pradeep.subedi@rutgers.edu
 */
#include <memory>

#include "DataSpacesReader.h"
#include "DataSpacesReader.tcc"
#include "adios2/helper/adiosCommMPI.h"
#include "adios2/helper/adiosFunctions.h" //CSVToVector
#include "adios2/toolkit/dataspaces/DSpacesConfig.h"
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

DataSpacesReader::DataSpacesReader(IO &io, const std::string &name, const Mode mode,
                                   helper::Comm comm)
: Engine("DataSpacesReader", io, name, mode, std::move(comm))
{

    f_Name = name;
    int ret = 0;
    latestStep = 0;
    nVars = 0;
    m_CurrentStep = -1;
    auto appID = m_IO.m_Parameters.find("AppID");
    if (appID != m_IO.m_Parameters.end())
    {
        m_data.appid = std::stoi(appID->second);
    }
    else
    {
        m_data.appid = 0;
    }
    auto latest = m_IO.m_Parameters.find("AlwaysProvideLatestTimestep");
    if (latest != m_IO.m_Parameters.end() && (latest->second == "yes" || latest->second == "true"))
    {
        m_ProvideLatest = true;
    }
    else
    {
        m_ProvideLatest = false;
    }
    MPI_Comm mpiComm = CommAsMPI(m_Comm);
    ret = adios_read_dataspaces_init(&mpiComm, &m_data);
    if (ret < 0)
    {
        fprintf(stderr, "DataSpaces did not initialize properly %d\n", ret);
    }
    else
    {
        m_IsOpen = true;
    }
}

DataSpacesReader::~DataSpacesReader()
{
    if (m_IsOpen)
    {
        DestructorClose(m_FailVerbose);
    }
    m_IsOpen = false;
}

StepStatus DataSpacesReader::BeginStep(StepMode mode, const float timeout_sec)
{
    // acquire lock, current step and n_vars in the Begin Step

    std::string ds_file_var, local_file_var;
    int elemsize, ndim, rank;
    int bcast_array[2] = {0, 0};
    uint64_t gdims[MAX_DS_NDIM], lb[MAX_DS_NDIM], ub[MAX_DS_NDIM];

    char *cstr, *meta_lk, *buffer;

    MPI_Comm_rank(m_data.mpi_comm, &rank);
    char *fstr = new char[f_Name.length() + 1];
    strcpy(fstr, f_Name.c_str());

    std::string lk_name = f_Name + std::to_string(m_CurrentStep + 1);
    meta_lk = new char[lk_name.length() + 1];
    strcpy(meta_lk, lk_name.c_str());

    MPI_Comm lock_comm = m_data.mpi_comm;
#ifndef HAVE_DSPACES2
    dspaces_lock_on_read(meta_lk, &lock_comm);
#endif /* HAVE_DSPACES2 */

    int nVars = 0;
    if (!m_ProvideLatest)
    {
        if (rank == 0)
        {
#ifdef HAVE_DSPACES2
            dspaces_client_t *client = get_client_handle();
            char meta_str[256];
            unsigned int metalen;
            snprintf(meta_str, sizeof(meta_str), "VARMETA@%s", fstr);
            int err = dspaces_get_meta(*client, meta_str, META_MODE_NEXT, m_CurrentStep,
                                       &bcast_array[1], (void **)&buffer, &metalen);
            bcast_array[0] = metalen;
#else
            buffer = dspaces_get_next_meta(m_CurrentStep, fstr, &bcast_array[0], &bcast_array[1]);
#endif /* HAVE_DSPACES2 */
        }
    }
    else
    {
        if (rank == 0)
        {
#ifdef HAVE_DSPACES2
            dspaces_client_t *client = get_client_handle();
            char meta_str[256];
            unsigned int metalen;
            snprintf(meta_str, sizeof(meta_str), "VARMETA@%s", fstr);
            int err = dspaces_get_meta(*client, meta_str, META_MODE_LAST, m_CurrentStep,
                                       &bcast_array[1], (void **)&buffer, &metalen);
            bcast_array[0] = metalen;
#else
            buffer = dspaces_get_latest_meta(m_CurrentStep, fstr, &bcast_array[0], &bcast_array[1]);
#endif /* HAVE_DSPACES2 */
        }
    }
    MPI_Bcast(bcast_array, 2, MPI_INT, 0, m_data.mpi_comm);
    int buf_len = bcast_array[0];
    int var_name_max_length = 128;
    nVars = (buf_len) /
            (2 * sizeof(int) + MAX_DS_NDIM * sizeof(uint64_t) + var_name_max_length * sizeof(char));
    m_CurrentStep = bcast_array[1];
    if (nVars == 0)
        return StepStatus::EndOfStream;
    if (rank != 0)
        buffer = (char *)malloc(buf_len);

    m_IO.RemoveAllVariables();

    MPI_Bcast(buffer, buf_len, MPI_CHAR, 0, m_data.mpi_comm);
    // now populate data from the buffer

    int *dim_meta, *elemSize_meta;
    uint64_t *gdim_meta;
    dim_meta = (int *)malloc(nVars * sizeof(int));
    elemSize_meta = (int *)malloc(nVars * sizeof(int));
    gdim_meta = (uint64_t *)malloc(MAX_DS_NDIM * nVars * sizeof(uint64_t));
    memset(gdim_meta, 0, MAX_DS_NDIM * nVars * sizeof(uint64_t));

    memcpy(dim_meta, buffer, nVars * sizeof(int));
    memcpy(elemSize_meta, &buffer[nVars * sizeof(int)], nVars * sizeof(int));
    memcpy(gdim_meta, &buffer[2 * nVars * sizeof(int)], MAX_DS_NDIM * nVars * sizeof(uint64_t));

    for (int var = 0; var < nVars; ++var)
    {

        int var_dim_size = dim_meta[var];
        int varType = elemSize_meta[var];

        std::string adiosName;
        char *val = (char *)(calloc(var_name_max_length, sizeof(char)));
        memcpy(val,
               &buffer[nVars * sizeof(int) + nVars * sizeof(int) +
                       nVars * MAX_DS_NDIM * sizeof(uint64_t) + var * var_name_max_length],
               var_name_max_length * sizeof(char));
        adiosName.assign(val);
        free(val);

        Dims shape;
        shape.resize(var_dim_size);
        if (var_dim_size > 0)
        {
            bool isOrderC = m_IO.m_ArrayOrder == adios2::ArrayOrdering::RowMajor;
            for (int i = 0; i < var_dim_size; i++)
            {
                if (isOrderC)
                {
                    shape[var_dim_size - i - 1] =
                        static_cast<int>(gdim_meta[var * MAX_DS_NDIM + i]);
                }
                else
                {
                    shape[i] = static_cast<int>(gdim_meta[var * MAX_DS_NDIM + i]);
                    ;
                }
            }
        }
        std::string adiosVarType;
        auto itType = ds_to_varType.find(varType);
        if (itType == ds_to_varType.end())
        {
            fprintf(stderr, "Wrong value in the serialized meta buffer for varType\n");
        }
        else
        {
            adiosVarType = itType->second;
        }
        if (adiosVarType == "int8_t")
            AddVar<int8_t>(m_IO, adiosName, shape);
        else if (adiosVarType == "uint8_t")
            AddVar<uint8_t>(m_IO, adiosName, shape);
        else if (adiosVarType == "int16_t")
            AddVar<int16_t>(m_IO, adiosName, shape);
        else if (adiosVarType == "uint16_t")
            AddVar<uint16_t>(m_IO, adiosName, shape);
        else if (adiosVarType == "int32_t")
            AddVar<int32_t>(m_IO, adiosName, shape);
        else if (adiosVarType == "uint32_t")
            AddVar<uint32_t>(m_IO, adiosName, shape);
        else if (adiosVarType == "int64_t")
            AddVar<int64_t>(m_IO, adiosName, shape);
        else if (adiosVarType == "uint64_t")
            AddVar<uint64_t>(m_IO, adiosName, shape);
        else if (adiosVarType == "float")
            AddVar<float>(m_IO, adiosName, shape);
        else if (adiosVarType == "double")
            AddVar<double>(m_IO, adiosName, shape);
        else if (adiosVarType == "long double")
            AddVar<long double>(m_IO, adiosName, shape);
        else if (adiosVarType == "float complex")
            AddVar<std::complex<float>>(m_IO, adiosName, shape);
        else if (adiosVarType == "double complex")
            AddVar<std::complex<double>>(m_IO, adiosName, shape);
        else
            AddVar<std::string>(m_IO, adiosName,
                                shape); // used string for last value
    }

    free(dim_meta);
    free(elemSize_meta);
    free(gdim_meta);
    free(buffer);

    return StepStatus::OK;
}

size_t DataSpacesReader::CurrentStep() const { return m_CurrentStep; }

void DataSpacesReader::EndStep()
{

    PerformGets();
    char *meta_lk;
    std::string lk_name = f_Name + std::to_string(m_CurrentStep);
    meta_lk = new char[lk_name.length() + 1];
    strcpy(meta_lk, lk_name.c_str());

    MPI_Comm lock_comm = m_data.mpi_comm;
#ifndef HAVE_DSPACES2
    dspaces_unlock_on_read(meta_lk, &lock_comm);
#endif /* HAVE_DSPACES2 */
}

void DataSpacesReader::DoClose(const int transportIndex)
{
    globals_adios_set_dataspaces_disconnected_from_reader();
}

void DataSpacesReader::Flush(const int transportIndex) {}

void DataSpacesReader::PerformGets()
{
    if (m_DeferredStack.size() > 0)
    {
#define declare_type(T)                                                                            \
    for (std::string variableName : m_DeferredStack)                                               \
    {                                                                                              \
        Variable<T> *var = m_IO.InquireVariable<T>(variableName);                                  \
        if (var != nullptr)                                                                        \
        {                                                                                          \
            ReadDsData(*var, var->GetData(), m_CurrentStep);                                       \
        }                                                                                          \
    }
        ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type
        m_DeferredStack.clear();
    }
}

#define declare_type(T)                                                                            \
    void DataSpacesReader::DoGetSync(Variable<T> &variable, T *data)                               \
    {                                                                                              \
        GetSyncCommon(variable, data);                                                             \
    }                                                                                              \
    void DataSpacesReader::DoGetDeferred(Variable<T> &variable, T *data)                           \
    {                                                                                              \
        GetDeferredCommon(variable, data);                                                         \
    }

ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

} // end namespace engine
} // end namespace core
} // end namespace adios2
