#ifndef ADIOS2_BINDINGS_CXX11_VARIABLE_DERIVED_H_
#define ADIOS2_BINDINGS_CXX11_VARIABLE_DERIVED_H_

#include "Operator.h"
#include "adios2/common/ADIOSTypes.h"

namespace adios2
{

/// \cond EXCLUDE_FROM_DOXYGEN
// forward declare
class IO; // friend
namespace core
{

class VariableDerived; // private implementation
}
/// \endcond

class VariableDerived
{
    friend class IO;

public:
    /**
     * Empty (default) constructor, use it as a placeholder for future
     * variables from IO:DefineVariableDerived<T> or IO:InquireVariableDerived<T>.
     * Can be used with STL containers.
     */
    VariableDerived() = default;

    /** Default, using RAII STL containers */
    ~VariableDerived() = default;

private:
    core::VariableDerived *m_VariableDerived = nullptr;

    VariableDerived(core::VariableDerived *variable);
};

} // end namespace adios2

#endif // ADIOS2_BINDINGS_CXX11_VARIABLE_DERIVED_H_
