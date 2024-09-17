/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * HDF5WriterP.h
 *
 *  Created on: March 20, 2017
 *      Author: Junmin
 */

#ifndef ADIOS2_ENGINE_HDF5_HDF5WRITERP_H__
#define ADIOS2_ENGINE_HDF5_HDF5WRITERP_H__

#include <hdf5.h>

#include "adios2/common/ADIOSConfig.h"
#include "adios2/core/Engine.h"
#include "adios2/core/IO.h"
#include "adios2/helper/adiosComm.h"
#include "adios2/toolkit/interop/hdf5/HDF5Common.h"

namespace adios2
{
namespace core
{
namespace engine
{

class HDF5WriterP : public Engine
{

public:
    /**
     * Constructor for HDF5 writer engine, writes in hdf5 format
     * @param name unique name given to the engine
     * @param accessMode
     * @param comm
     * @param method
     */
    HDF5WriterP(IO &io, const std::string &name, const Mode mode, helper::Comm comm);

    ~HDF5WriterP();

    StepStatus BeginStep(StepMode mode, const float timeoutSeconds = -1.0) final;
    void EndStep() final;

    void PerformPuts() final;

private:
    interop::HDF5Common m_H5File;
    bool m_Flushed = false;

    void Init();

#define declare_type(T)                                                                            \
    void DoPutSync(Variable<T> &variable, const T *values) final;                                  \
    void DoPutDeferred(Variable<T> &variable, const T *values) final;
    ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

    template <class T>
    void DoPutSyncCommon(Variable<T> &variable, const T *values);

    void DoClose(const int transportIndex = -1) final;

    /**
     * Called if destructor is called on an open engine.  Should warn or take
     * any non-complex measure that might help recover.
     */
    void DestructorClose(bool Verbose) noexcept final{};

    void Flush(const int transportIndex = -1) final;
};

} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_ENGINE_HDF5_HDF5WRITERP_H__ */
