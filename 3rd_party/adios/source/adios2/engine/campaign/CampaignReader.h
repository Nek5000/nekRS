/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * CampaignReader.h
 * An empty skeleton engine from which any engine can be built
 *
 *  Created on: May 15, 2023
 *      Author: Norbert Podhorszki pnorbert@ornl.gov
 */

#ifndef ADIOS2_ENGINE_CAMPAIGNREADER_H_
#define ADIOS2_ENGINE_CAMPAIGNREADER_H_

#include "CampaignData.h"
#include "adios2/common/ADIOSConfig.h"
#include "adios2/core/ADIOS.h"
#include "adios2/core/Engine.h"
#include "adios2/helper/adiosComm.h"
#include "adios2/helper/adiosFunctions.h"

#include <sqlite3.h>

namespace adios2
{
namespace core
{
namespace engine
{

class CampaignReader : public Engine
{
public:
    /**
     * Constructor for single BP capsule engine, writes in BP format into a
     * single
     * heap capsule
     * @param name unique name given to the engine
     * @param accessMode
     * @param comm
     * @param method
     * @param hostLanguage
     */
    CampaignReader(IO &adios, const std::string &name, const Mode mode, helper::Comm comm);

    ~CampaignReader();
    StepStatus BeginStep(StepMode mode = StepMode::Read, const float timeoutSeconds = -1.0) final;
    void PerformGets() final;
    size_t CurrentStep() const final;
    void EndStep() final;

    MinVarInfo *MinBlocksInfo(const VariableBase &, const size_t Step) const;
    bool VarShape(const VariableBase &Var, const size_t Step, Dims &Shape) const;
    bool VariableMinMax(const VariableBase &, const size_t Step, MinMaxStruct &MinMax);

private:
    UserOptions::Campaign m_Options;
    int m_ReaderRank; // my rank in the readers' comm

    int m_CurrentStep = 0;

    // EndStep must call PerformGets if necessary
    bool m_NeedPerformGets = false;

    std::vector<adios2::core::IO *> m_IOs;
    std::vector<adios2::core::Engine *> m_Engines;

    struct VarInternalInfo
    {
        void *originalVar; // Variable<T> in the actual IO
        size_t ioIdx;      // actual IO in m_IOs
        size_t engineIdx;  // actual engine in m_Engines
        VarInternalInfo(void *p, size_t i, size_t e) : originalVar(p), ioIdx(i), engineIdx(e) {}
    };
    std::unordered_map<std::string, VarInternalInfo> m_VarInternalInfo;

    void Init() final; ///< called from constructor, gets the selected Skeleton
                       /// transport method from settings
    void ReadConfig(std::string path);
    void InitParameters() final;
    void InitTransports() final;

#define declare_type(T)                                                                            \
    void DoGetSync(Variable<T> &, T *) final;                                                      \
    void DoGetDeferred(Variable<T> &, T *) final;
    ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

    void DoClose(const int transportIndex = -1);

#define declare_type(T)                                                                            \
    std::map<size_t, std::vector<typename Variable<T>::BPInfo>> DoAllStepsBlocksInfo(              \
        const Variable<T> &variable) const final;                                                  \
                                                                                                   \
    std::vector<std::vector<typename Variable<T>::BPInfo>> DoAllRelativeStepsBlocksInfo(           \
        const Variable<T> &) const final;                                                          \
                                                                                                   \
    std::vector<typename Variable<T>::BPInfo> DoBlocksInfo(const Variable<T> &variable,            \
                                                           const size_t step) const final;

    ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

    /**
     * Called if destructor is called on an open engine.  Should warn or take
     * any non-complex measure that might help recover.
     */
    void DestructorClose(bool Verbose) noexcept final;

    /**
     * Create a new variable with name `name` in `io`
     * based on an existing variable.
     */
    template <class T>
    Variable<T> DuplicateVariable(Variable<T> *variable, IO &io, std::string &name,
                                  VarInternalInfo &vii);

    /**
     * Create a new attribute with name `name` in `io`
     * based on an existing attribute.
     */
    template <class T>
    Attribute<T> DuplicateAttribute(Attribute<T> *attribute, IO &io, std::string &name);

    /**
     * Create a new variable with name `name` in `io`
     * based on an existing variable.
     */
    template <class T>
    std::pair<Variable<T> *, Engine *> TranslateToActualVariable(Variable<T> &variable);

    sqlite3 *m_DB;
    CampaignData m_CampaignData;
};

} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_ENGINE_CAMPAIGNREADER_H_ */
