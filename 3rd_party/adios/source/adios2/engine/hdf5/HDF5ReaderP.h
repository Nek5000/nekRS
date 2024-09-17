/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * HDF5ReaderP.h
 *
 *  Created on: March 20, 2017
 *      Author: Junmin
 */

#ifndef ADIOS2_ENGINE_HDF5_HDF5READERP_H_
#define ADIOS2_ENGINE_HDF5_HDF5READERP_H_

#include "adios2/core/Engine.h"
#include "adios2/core/IO.h"
#include "adios2/helper/adiosComm.h"
#include "adios2/toolkit/interop/hdf5/HDF5Common.h"

#include <map>
#include <vector>

namespace adios2
{
namespace core
{
namespace engine
{

class HDF5ReaderP : public Engine
{

public:
    /**
     * Constructor for single HDF5 reader engine, reads from HDF5 format
     * @param name unique name given to the engine
     * @param accessMode
     * @param comm
     * @param method
     */
    HDF5ReaderP(IO &adios, const std::string &name, const Mode openMode, helper::Comm comm);

    ~HDF5ReaderP();

    bool IsValid();

    StepStatus BeginStep(StepMode mode, const float timeoutSeconds = -1.0) final;
    size_t CurrentStep() const final;
    void EndStep() final;

    void PerformGets() final;

private:
    interop::HDF5Common m_H5File;
    void Init() final;

    bool m_InStreamMode = false; // default is not streaming, i.e. set var timestep range
                                 // unsigned int m_StreamAt = 0; // stream step counter
    size_t m_StreamAt = 0;
#define declare_type(T)                                                                            \
    void DoGetSync(Variable<T> &, T *) final;                                                      \
    void DoGetDeferred(Variable<T> &, T *) final;                                                  \
    std::map<size_t, std::vector<typename Variable<T>::BPInfo>> DoAllStepsBlocksInfo(              \
        const Variable<T> &variable) const final;                                                  \
                                                                                                   \
    std::vector<typename Variable<T>::BPInfo> DoBlocksInfo(const Variable<T> &variable,            \
                                                           const size_t step) const final;
    ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

    void DoClose(const int transportIndex = -1) final;

    template <class T>
    size_t ReadDataset(hid_t dataSetId, hid_t h5Type, Variable<T> &variable, T *values);

    template <class T>
    void GetSyncCommon(Variable<T> &variable, T *data);

    template <class T>
    void GetDeferredCommon(Variable<T> &variable, T *data);

    template <class T>
    std::map<size_t, std::vector<typename Variable<T>::BPInfo>>
    GetAllStepsBlocksInfo(const Variable<T> &variable) const;

    template <class T>
    std::vector<typename Variable<T>::BPInfo> GetBlocksInfo(const Variable<T> &variable,
                                                            const size_t step) const;

    template <class T>
    std::vector<typename core::Variable<T>::BPInfo>
    BlocksInfoCommon(const core::Variable<T> &variable) const;

    template <class T>
    void UseHDFRead(Variable<T> &variable, T *values, hid_t h5Type);

    std::vector<std::string> m_DeferredStack;

    size_t DoSteps() const final;
};

} // end namespace engine
} // end namespace core
} // end namespace adios2
#endif /* ADIOS2_ENGINE_HDF5_HDF5READERP_H_ */
