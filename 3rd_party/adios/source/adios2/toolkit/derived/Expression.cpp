#ifndef ADIOS2_DERIVED_Expression_CPP_
#define ADIOS2_DERIVED_Expression_CPP_

#include "Expression.h"
#include "parser/ASTDriver.h"

namespace adios2
{
namespace detail
{
// helper function
adios2::detail::ExpressionOperator convert_op(std::string opname)
{
    adios2::detail::ExpressionOperator op;
    try
    {
        op = adios2::detail::get_op(opname);
    }
    catch (std::out_of_range &e)
    {
        (void)e; // use e
        helper::Throw<std::invalid_argument>("Derived", "ExprHelper", "get_op",
                                             "Parser cannot recognize operator '" + opname + "'.");
    }
    return op;
};

adios2::derived::ExpressionTree ASTNode_to_ExpressionTree(adios2::detail::ASTNode *node)
{
    adios2::derived::ExpressionTree exprTree_node(convert_op(node->get_opname()));
    for (adios2::detail::ASTNode *e : node->get_subexprs())
    {
        switch (convert_op(e->get_opname()))
        {
        case adios2::detail::ExpressionOperator::OP_ALIAS: // add variable given by alias
            // add an index operation in the chain if the variable contains indeces
            /*if (e->lookup_var_indices(e->alias) != "")
            {
                ExpressionTree index_expr(adios2::detail::ExpressionOperator::OP_INDEX);
                index_expr.set_indeces(e->lookup_var_indices(e->alias));
                index_expr.add_child(e->lookup_var_path(e->alias));
                expTree_node->add_child(expr);
            }*/
            exprTree_node.add_child(e->get_varname());
            break;
        case adios2::detail::ExpressionOperator::OP_PATH: // add variable name
            exprTree_node.add_child(e->get_varname());
            break;
        case adios2::detail::ExpressionOperator::OP_NUM: // set the base value for the operation
            exprTree_node.set_base(e->get_value());
            break;
        default: // if the children nodes are other expressions, convert them to expressions
            auto temp_node = ASTNode_to_ExpressionTree(e);
            // move from a binary to a multinary tree if the child has the same operation
            if (convert_op(e->get_opname()) == convert_op(node->get_opname()) &&
                adios2::detail::op_property.at(convert_op(e->get_opname())).is_associative)
            {
                // concatenate exprTree with temp_node
                for (std::tuple<adios2::derived::ExpressionTree, std::string, bool> childTree :
                     temp_node.sub_exprs)
                {
                    if (std::get<2>(childTree) == true)
                        exprTree_node.add_child(std::get<0>(childTree));
                    else
                        exprTree_node.add_child(std::get<1>(childTree));
                }
            }
            else
            {
                exprTree_node.add_child(temp_node);
            }
        }
    }
    return exprTree_node;
}
}

namespace derived
{

Expression::Expression(std::string string_exp)
: m_Shape({0}), m_Start({0}), m_Count({0}), ExprString(string_exp)
{
    adios2::detail::ASTDriver drv(string_exp);
    m_Expr = adios2::detail::ASTNode_to_ExpressionTree(drv.getAST());
}

std::vector<std::string> Expression::VariableNameList() { return m_Expr.VariableNameList(); }

Dims Expression::GetShape() { return m_Shape; }

Dims Expression::GetStart() { return m_Start; }

Dims Expression::GetCount() { return m_Count; }

void Expression::SetDims(std::map<std::string, std::tuple<Dims, Dims, Dims>> NameToDims)
{
    std::map<std::string, Dims> NameToCount, NameToStart, NameToShape;
    for (const auto &it : NameToDims)
    {
        NameToStart[it.first] = std::get<0>(it.second);
        NameToCount[it.first] = std::get<1>(it.second);
        NameToShape[it.first] = std::get<2>(it.second);
    }
    m_Count = m_Expr.GetDims(NameToCount);
    m_Start = m_Expr.GetDims(NameToStart);
    m_Shape = m_Expr.GetDims(NameToShape);
}

std::vector<DerivedData>
Expression::ApplyExpression(DataType type, size_t numBlocks,
                            std::map<std::string, std::vector<DerivedData>> nameToData)
{
    return m_Expr.ApplyExpression(type, numBlocks, nameToData);
}

void ExpressionTree::set_base(double c) { detail.constant = c; }

void ExpressionTree::set_indeces(std::vector<std::tuple<size_t, size_t, size_t>> index_list)
{
    detail.indices = index_list;
}

void ExpressionTree::add_child(ExpressionTree exp) { sub_exprs.push_back({exp, "", true}); }

void ExpressionTree::add_child(std::string var)
{
    sub_exprs.push_back({ExpressionTree(), var, false});
}

std::vector<std::string> ExpressionTree::VariableNameList()
{
    std::vector<std::string> var_list;
    for (auto subexp : sub_exprs)
    {
        // if the sub_expression is a leaf
        if (!std::get<2>(subexp))
        {
            var_list.push_back(std::get<1>(subexp));
        }
        else
        {
            auto subexpr_list = std::get<0>(subexp).VariableNameList();
            var_list.insert(var_list.end(), subexpr_list.begin(), subexpr_list.end());
        }
    }
    return var_list;
}

void ExpressionTree::print()
{
    std::cout << "Print Expression:" << std::endl;
    std::cout << "\toperation: " << get_op_name(detail.operation) << std::endl;
    std::cout << "\tconstant: " << detail.constant << std::endl;
    std::cout << "\tchildren: " << sub_exprs.size() << std::endl;

    for (std::tuple<ExpressionTree, std::string, bool> t : sub_exprs)
    {
        if (std::get<2>(t) == true)
        {
            std::get<0>(t).print();
        }
        else
        {
            std::cout << "string: " << std::get<1>(t) << std::endl;
        }
    }
}

Dims ExpressionTree::GetDims(std::map<std::string, Dims> NameToDims)
{
    std::vector<Dims> exprDims;
    for (auto subexp : sub_exprs)
    {
        // if the sub_expression is a leaf, we get the shape from the input std::map
        if (!std::get<2>(subexp))
        {
            exprDims.push_back(NameToDims[std::get<1>(subexp)]);
        }
        else
        {
            exprDims.push_back(std::get<0>(subexp).GetDims(NameToDims));
        }
    }
    // get the output dimensions after applying the operator
    auto op_fct = OpFunctions.at(detail.operation);
    Dims opDims = op_fct.DimsFct(exprDims);
    return opDims;
}

std::vector<DerivedData>
ExpressionTree::ApplyExpression(DataType type, size_t numBlocks,
                                std::map<std::string, std::vector<DerivedData>> nameToData)
{
    // create operands for the computation function
    // exprData[0] = list of void* data for block 0 for each variable
    std::vector<std::vector<DerivedData>> exprData(numBlocks);
    std::vector<bool> dealocate;
    for (auto subexp : sub_exprs)
    {
        if (!std::get<2>(subexp))
        {
            // do not dealocate leafs (this is user data)
            dealocate.push_back(false);
            for (size_t blk = 0; blk < numBlocks; blk++)
            {
                exprData[blk].push_back(nameToData[std::get<1>(subexp)][blk]);
            }
        }
        else
        {
            dealocate.push_back(true);
            auto subexpData = std::get<0>(subexp).ApplyExpression(type, numBlocks, nameToData);
            for (size_t blk = 0; blk < numBlocks; blk++)
            {
                exprData[blk].push_back(subexpData[blk]);
            }
        }
    }
    // apply the computation operator on all blocks
    std::vector<DerivedData> outputData(numBlocks);
    auto op_fct = OpFunctions.at(detail.operation);
    for (size_t blk = 0; blk < numBlocks; blk++)
    {
        outputData[blk] = op_fct.ComputeFct(exprData[blk], type);
    }
    // deallocate intermediate data after computing the operation
    for (size_t blk = 0; blk < numBlocks; blk++)
    {
        for (size_t i = 0; i < exprData[blk].size(); i++)
        {
            if (dealocate[i] == false)
                continue;
            free(exprData[blk][i].Data);
        }
    }
    return outputData;
}

}
}
#endif
