/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * SstWriter.h
 *
 *  Created on: Aug 17, 2017
 *      Author: Greg Eisenhauer
 */

#ifndef ADIOS2_ENGINE_SST_SST_WRITER_H_
#define ADIOS2_ENGINE_SST_SST_WRITER_H_

#include "adios2/common/ADIOSConfig.h"
#include "adios2/core/Engine.h"
#include "adios2/helper/adiosComm.h"
#include "adios2/toolkit/format/bp/bp3/BP3Serializer.h"
#include "adios2/toolkit/format/bp5/BP5Serializer.h"
#include "adios2/toolkit/sst/sst.h"

#include <memory>

namespace adios2
{
namespace core
{
namespace engine
{

class SstWriter : public Engine
{

public:
    SstWriter(IO &io, const std::string &name, const Mode mode, helper::Comm comm);

    virtual ~SstWriter();

    StepStatus BeginStep(StepMode mode, const float timeoutSeconds = -1.0) final;
    size_t CurrentStep() const final;
    void PerformPuts() final;
    void EndStep() final;
    void Flush(const int transportIndex = -1) final;
    void NotifyEngineAttribute(std::string name, DataType type) noexcept;
    void NotifyEngineAttribute(std::string name, AttributeBase *Attr, void *data) noexcept;

    /**
     * Called if destructor is called on an open engine.  Should warn or take
     * any non-complex measure that might help recover.
     */
    void DestructorClose(bool Verbose) noexcept final;

private:
    void Init(); ///< calls InitCapsules and InitTransports based on Method,
                 /// called from constructor

#define declare_type(T)                                                                            \
    void DoPutSync(Variable<T> &variable, const T *values) final;                                  \
    void DoPutDeferred(Variable<T> &, const T *) final;
    ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

    void DoPutStructSync(VariableStruct &, const void *) final;
    void DoPutStructDeferred(VariableStruct &, const void *) final;

    void PutStructCommon(VariableBase &, const void *);

    template <class T>
    void PutSyncCommon(Variable<T> &variable, const T *values);

    template <class T>
    void PutDeferredCommon(Variable<T> &variable, const T *values);

    struct BP3DataBlock
    {
        _SstData data;
        _SstData metadata;
        format::BP3Serializer *serializer;
    };
    struct BP5DataBlock
    {
        _SstData data;
        _SstData metadata;
        _SstData attribute_data;
        SstMetaMetaList MetaMetaBlocks;
        format::BP5Serializer::TimestepInfo *TSInfo;
    };
    std::unique_ptr<format::BP3Serializer> m_BP3Serializer;

    std::unique_ptr<format::BP5Serializer> m_BP5Serializer;

    SstStream m_Output;
    long m_WriterStep = -1;
    bool m_DefinitionsNotified = false;
    bool m_MarshalAttributesNecessary = true; // first time through, marshal
    struct _SstParams Params;

    void MarshalAttributes();
    void DoClose(const int transportIndex = -1) final;
};

} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_ENGINE_SST_SST_WRITER_H_ */
