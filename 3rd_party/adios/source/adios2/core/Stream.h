/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Stream.h : higher-level fstream-like API
 *
 *  Created on: Jan 5, 2018
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_CORE_STREAM_H_
#define ADIOS2_CORE_STREAM_H_

#include <memory> //std::shared_ptr

#include "adios2/common/ADIOSMacros.h"
#include "adios2/common/ADIOSTypes.h"
#include "adios2/core/ADIOS.h"
#include "adios2/core/Engine.h"
#include "adios2/core/IO.h"
#include "adios2/helper/adiosComm.h"

namespace adios2
{
namespace core
{

/**
 * Private implementation of high-level APIs
 */
class Stream
{

public:
    /** internal ADIOS object owned by Stream */
    std::shared_ptr<ADIOS> m_ADIOS;

    /**
     * reference to IO with name m_Name inside m_ADIOS, exposed as public
     * to use certain functions
     */
    IO *m_IO = nullptr;

    /**
     * reference to Engine with name m_Name inside m_ADIOS, exposed as public
     * to use certain functions directly
     */
    Engine *m_Engine = nullptr;

    Stream(const std::string &name, const Mode mode, helper::Comm comm,
           const std::string engineType, const std::string hostLanguage);

    Stream(const std::string &name, const Mode mode, const std::string engineType,
           const std::string hostLanguage);

    Stream(const std::string &name, const Mode mode, helper::Comm comm,
           const std::string configFile, const std::string ioInConfigFile,
           const std::string hostLanguage);

    Stream(const std::string &name, const Mode mode, const std::string configFile,
           const std::string ioInConfigFile, const std::string hostLanguage);

    ~Stream() = default;

    template <class T>
    void WriteAttribute(const std::string &name, const T &value, const std::string &variableName,
                        const std::string separator, const bool nextStep);
    template <class T>
    void WriteAttribute(const std::string &name, const T *array, const size_t elements,
                        const std::string &variableName, const std::string separator,
                        const bool nextStep);

    template <class T>
    void Write(const std::string &name, const T *values, const Dims &shape, const Dims &start,
               const Dims &count, const vParams &operations = vParams(),
               const bool endStep = false);

    template <class T>
    void Write(const std::string &name, const T &value, const bool isLocalValue,
               const bool endStep);

    bool GetStep();

    template <class T>
    void Read(const std::string &name, T *values, const size_t blockID);

    template <class T>
    void Read(const std::string &name, T *values, const Box<size_t> &step, const size_t blockID);

    template <class T>
    void Read(const std::string &name, T *values, const Box<Dims> &selection, const size_t blockID);

    template <class T>
    void Read(const std::string &name, T *values, const Box<Dims> &selection,
              const Box<size_t> &stepSelection, const size_t blockID);

    template <class T>
    std::vector<T> Read(const std::string &name, const size_t blockID);

    template <class T>
    std::vector<T> Read(const std::string &name, const Box<size_t> &stepsSelection,
                        const size_t blockID);

    template <class T>
    std::vector<T> Read(const std::string &name, const Box<Dims> &selection, const size_t blockID);

    template <class T>
    std::vector<T> Read(const std::string &name, const Box<Dims> &selection,
                        const Box<size_t> &stepsSelection, const size_t blockID);

    template <class T>
    void ReadAttribute(const std::string &name, T *data, const std::string &variableName,
                       const std::string separator);

    void EndStep();

    void Close();

    size_t CurrentStep() const;

    size_t Steps() const;

private:
    /** Stream, IO and Engine names  */
    std::string m_Name;

    /** mode to open engine at first read/write */
    Mode m_Mode;

    /** Sets engine type to be opened at first read/write */
    std::string m_EngineType;

    /** internal flag to check if getstep was called */
    bool m_FirstStep = true;

    bool m_StepStatus = false;

    template <class T>
    std::vector<T> GetCommon(Variable<T> &variable);

    template <class T>
    void GetPCommon(Variable<T> &variable, T *values);

    template <class T>
    void CheckPCommon(const std::string &name, const T *values) const;

    void CheckOpen();

    template <class T>
    void SetBlockSelectionCommon(Variable<T> &variable, const size_t blockID);
};

} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_CORE_STREAM_H_ */
