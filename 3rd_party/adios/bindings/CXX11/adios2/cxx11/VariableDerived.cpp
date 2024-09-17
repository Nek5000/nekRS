#include "VariableDerived.h"

#include "adios2/core/VariableDerived.h"

namespace adios2
{
VariableDerived::VariableDerived(core::VariableDerived *variable) : m_VariableDerived(variable) {}
} // end namespace adios2
