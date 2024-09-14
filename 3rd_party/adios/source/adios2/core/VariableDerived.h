#ifndef ADIOS2_CORE_VARIABLE_DERIVED_H_
#define ADIOS2_CORE_VARIABLE_DERIVED_H_

#include "adios2/common/ADIOSTypes.h"
#include "adios2/core/VariableBase.h"
#include "adios2/helper/adiosType.h"
#include "adios2/toolkit/derived/Expression.h"

namespace adios2
{
namespace core
{

/**
 * @param Base (parent) class for template derived (child) class Variable.
 */
class VariableDerived : public VariableBase
{
    DerivedVarType m_DerivedType;

public:
    adios2::derived::Expression m_Expr;
    VariableDerived(const std::string &name, adios2::derived::Expression expr,
                    const DataType exprType, const bool isConstant, const DerivedVarType varType);
    ~VariableDerived() = default;

    DerivedVarType GetDerivedType();
    std::vector<std::string> VariableNameList();
    void UpdateExprDim(std::map<std::string, std::tuple<Dims, Dims, Dims>> NameToDims);

    std::vector<void *>
    ApplyExpression(std::map<std::string, std::vector<void *>> NameToData,
                    std::map<std::string, std::tuple<Dims, Dims, Dims>> NameToDims);
    std::vector<std::tuple<void *, Dims, Dims>>
    ApplyExpression(std::map<std::string, std::unique_ptr<MinVarInfo>> &mvi);
};

} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_CORE_VARIABLE_DERIVED_H_ */
