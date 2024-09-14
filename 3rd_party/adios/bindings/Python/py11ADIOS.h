/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * ADIOSPy.h  python binding to ADIOS class
 *
 *  Created on: Mar 13, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_BINDINGS_PYTHON_ADIOS_H_
#define ADIOS2_BINDINGS_PYTHON_ADIOS_H_

#include "py11IO.h"
#include "py11Operator.h"

#include <memory> //std::shared_ptr
#include <string>

#include "adios2/core/ADIOS.h"

namespace adios2
{
namespace py11
{

class ADIOS
{

public:
#if ADIOS2_USE_MPI
    ADIOS(const std::string &configFile, MPI4PY_Comm comm);
    ADIOS(MPI4PY_Comm comm);
#endif
    ADIOS(const std::string &configFile);
    ADIOS();
    ~ADIOS() = default;

    /** object inspection true: valid object, false: invalid object */
    explicit operator bool() const noexcept;

    IO DeclareIO(const std::string name);
    IO AtIO(const std::string name);

    Operator DefineOperator(const std::string name, const std::string type,
                            const Params &parameters = Params());

    Operator InquireOperator(const std::string name);

    bool RemoveIO(const std::string name);

    void RemoveAllIOs();

    void FlushAll();

private:
    std::shared_ptr<adios2::core::ADIOS> m_ADIOS;

    void CheckPointer(const std::string hint);
};

} // end namespace py11
} // end namespace adios2

#endif /* ADIOS2_BINDINGS_PYTHON_ADIOS_H_ */
