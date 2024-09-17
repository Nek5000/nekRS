/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * HDF5CommonP.h
 *
 *  Created on: March 20, 2017
 *      Author: Junmin
 */

#ifndef ADIOS2_TOOLKIT_INTEROP_HDF5_HDF5COMMON_H_
#define ADIOS2_TOOLKIT_INTEROP_HDF5_HDF5COMMON_H_

#include <hdf5.h>

#include <string>

#include "adios2/common/ADIOSMacros.h"
#include "adios2/common/ADIOSTypes.h"
#include "adios2/core/IO.h" // for CreateVar
#include "adios2/core/Variable.h"
#include "adios2/helper/adiosComm.h"

#include <stdexcept> // for Intel Compiler

namespace adios2
{
namespace interop
{

typedef enum
{
    E_H5_DATASET = 0,
    E_H5_DATATYPE = 1,
    E_H5_GROUP = 2,
    E_H5_SPACE = 3,
    E_H5_ATTRIBUTE = 4
} ADIOS_ENUM_H5;

class HDF5DatasetGuard
{
public:
    HDF5DatasetGuard(std::vector<hid_t> &chain) : m_Chain(chain) {}

    ~HDF5DatasetGuard()
    {
        if (m_Chain.size() == 0)
        {
            return;
        }

        for (size_t i = 0; i < m_Chain.size() - 1; i++)
        {
            H5Gclose(m_Chain[i]);
        }
        H5Dclose(m_Chain.back());
    }

private:
    std::vector<hid_t> m_Chain;
};

class HDF5TypeGuard
{
public:
    HDF5TypeGuard(hid_t key, ADIOS_ENUM_H5 type)
    {
        m_Key = key;
        m_Type = type;
        if (key < 0)
        {
            helper::Throw<std::ios_base::failure>("Toolkit",
                                                  "interop::hdf5::HDF5Common::HDF5TypeGuard",
                                                  "HDF5TypeGuard", "HDF5 failure detected");
        }
    }

    ~HDF5TypeGuard()
    {
        if (m_Type == E_H5_DATASET)
        {
            H5Dclose(m_Key);
        }
        else if (m_Type == E_H5_GROUP)
        {
            H5Gclose(m_Key);
        }
        else if (m_Type == E_H5_SPACE)
        {
            H5Sclose(m_Key);
        }
        else if (m_Type == E_H5_DATATYPE)
        {
            H5Tclose(m_Key);
        }
        else if (m_Type == E_H5_ATTRIBUTE)
        {
            H5Aclose(m_Key);
        }
        else
        {
            printf(" UNABLE to close \n");
        }
    }

private:
    ADIOS_ENUM_H5 m_Type;
    hid_t m_Key;
};

class HDF5Common
{

public:
    /**
     * Unique constructor for HDF5 file
     */
    HDF5Common();
    ~HDF5Common();

    static const std::string ATTRNAME_NUM_STEPS;
    static const std::string ATTRNAME_GIVEN_ADIOSNAME;
    static const std::string PREFIX_BLOCKINFO;
    static const std::string PREFIX_STAT;

    static const std::string PARAMETER_COLLECTIVE;
    static const std::string PARAMETER_CHUNK_FLAG;
    static const std::string PARAMETER_CHUNK_VARS;
    static const std::string PARAMETER_HAS_IDLE_WRITER_RANK;

    void ParseParameters(core::IO &io);
    void Init(const std::string &name, helper::Comm const &comm, bool toWrite);
    void Append(const std::string &name, helper::Comm const &comm);

    template <class T>
    void Write(core::Variable<T> &variable, const T *values);

    /*
     * This function will define a non string variable to HDF5
     * note that define a dataset in HDF5 means allocate space and place an
     * entry to the HDF5 file. (By default we define variable when a PUT is
     * called from adios client)
     */
    template <class T>
    void DefineDataset(core::Variable<T> &variable);

    void CreateDataset(const std::string &varName, hid_t h5Type, hid_t filespaceID,
                       std::vector<hid_t> &chain);
    bool OpenDataset(const std::string &varName, std::vector<hid_t> &chain);
    void RemoveEmptyDataset(const std::string &varName);
    void StoreADIOSName(const std::string adiosName, hid_t dsetID);
    void ReadADIOSName(hid_t dsetID, std::string &adiosName);

    void LocateAttrParent(const std::string &attrName, std::vector<std::string> &list,
                          std::vector<hid_t> &datasetChain);

    void Close();
    void Advance();

    /*
     * This function will browse all (non-string) variables in io and define
     * them in HDF5 This is a back up mode compare to the default behaviour that
     * defines a variable in HDF5 when the PUT is called on that variable. This
     * function is expected to be called by BeginStep(), for collectiveness,
     * required by HDF5
     */
    void CreateVarsFromIO(core::IO &io);
    void CleanUpNullVars(core::IO &io);

    void WriteAttrFromIO(core::IO &io);
    void ReadAttrToIO(core::IO &io);

    void ReadNativeAttrToIO(core::IO &io, hid_t datasetId, std::string const &pathFromRoot);

    void SetAdiosStep(size_t ts);

    size_t GetNumAdiosSteps();
    size_t GetAdiosStep() const;
    void WriteAdiosSteps();

    void ReadVariables(size_t ts, core::IO &io);
    void FindVarsFromH5(core::IO &io, hid_t gid, const char *name, const char *heritage, size_t ts);
    void ReadAllVariables(core::IO &io);

    void ReadStringScalarDataset(hid_t dataSetId, std::string &result);
    hid_t GetTypeStringScalar(const std::string &input);
    void CreateVar(core::IO &io, hid_t h5Type, std::string const &name, size_t ts);

    template <class T>
    void AddVar(core::IO &io, std::string const &name, hid_t datasetId, size_t ts);

    // adios only allows a scalar string var
    void AddSingleString(core::IO &io, std::string const &name, hid_t datasetId, size_t ts);

    // decompose array string vars
    void AddVarString(core::IO &io, std::string const &name, hid_t datasetId, size_t ts);

    template <class T>
    void AddNonStringAttribute(core::IO &io, std::string const &attrName, hid_t attrId,
                               hid_t h5Type, hsize_t arraySize);

    static void StaticGetAdiosStepString(std::string &adiosStepName, size_t ts);

    hid_t m_PropertyListId = -1;
    hid_t m_PropertyTxfID = -1;
    hid_t m_FileId = -1;
    hid_t m_GroupId = -1;

    hid_t m_DefH5TypeComplexLongDouble;
    hid_t m_DefH5TypeComplexDouble;
    hid_t m_DefH5TypeComplexFloat;
    hid_t m_DefH5TypeBlockStat;

    size_t m_CurrentAdiosStep = 0;

    void CheckWriteGroup();

    template <class T>
    hid_t GetHDF5Type(); // should this be public?

    template <class T>
    T GetADIOSType(hid_t);

    bool m_IsGeneratedByAdios = false;

    struct MPI_API
    {
        bool (*init)(helper::Comm const &comm, hid_t id, int *rank, int *size);
        herr_t (*set_dxpl_mpio)(hid_t dxpl_id, H5FD_mpio_xfer_t xfer_mode);
    };

private:
    void ReadInStringAttr(core::IO &io, const std::string &attrName, hid_t attrId, hid_t h5Type,
                          hid_t sid);
    void ReadInNonStringAttr(core::IO &io, const std::string &attrName, hid_t attrId, hid_t h5Type,
                             hid_t sid);
    void WriteStringAttr(core::IO &io, core::Attribute<std::string> *adiosAttr,
                         const std::string &attrName, hid_t parentID);

    template <class T>
    void WriteNonStringAttr(core::IO &io, core::Attribute<T> *adiosAttr, hid_t parentID,
                            const char *h5Name);

    template <class T>
    void GetHDF5SpaceSpec(const core::Variable<T> &variable, std::vector<hsize_t> &,
                          std::vector<hsize_t> &, std::vector<hsize_t> &);

    void CheckVariableOperations(const core::VariableBase &variable) const;

    bool m_WriteMode = false;

    size_t m_NumAdiosSteps = 0;

    MPI_API const *m_MPI = nullptr;
    int m_CommRank = 0;
    int m_CommSize = 1;

    template <class T>
    void AddBlockInfo(const core::Variable<T> &varaible, hid_t parentId);
    template <class T>
    void AddStats(const core::Variable<T> &variable, hid_t parentId, std::vector<T> &stats);

    hid_t m_TimeStepH5T = H5T_NATIVE_ULLONG;
    ;
    hid_t m_ChunkPID;
    size_t m_ChunkDim;
    std::set<std::string> m_ChunkVarNames;
    bool m_OrderByC = true; // C or fortran

    // Some write rank can be idle. This causes conflict with HDF5 collective
    // requirement in functions Guard this by load vars in beginStep
    bool m_IdleWriterOn = false;
};

} // end namespace interop
} // end namespace adios

#endif /* ADIOS2_TOOLKIT_INTEROP_HDF5_HDF5COMMON_H_ */
