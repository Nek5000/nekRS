#include "ASTNode.h"

namespace adios2
{
namespace detail
{

ASTNode::ASTNode() {}

ASTNode::ASTNode(std::string op) { opname = op; }

ASTNode::ASTNode(std::string op, size_t numsubexprs)
{
    opname = op;
    sub_exprs.resize(numsubexprs);
}

ASTNode::ASTNode(std::string op, std::string a)
{
    opname = op;
    alias = a;
}

ASTNode::ASTNode(std::string op, std::vector<std::tuple<int, int, int>> i)
{
    opname = op;
    indices = i;
}

ASTNode::~ASTNode()
{
    for (ASTNode *sub_expr : sub_exprs)
    {
        delete sub_expr;
    }
    sub_exprs.clear();
}

void ASTNode::set_num_subexprs(size_t n) { sub_exprs.resize(n); }

void ASTNode::pushback_subexpr(ASTNode *subexpr) { sub_exprs.push_back(subexpr); }

void ASTNode::insert_subexpr_n(ASTNode *subexpr, size_t index) { sub_exprs[index] = subexpr; }

std::string ASTNode::printpretty(std::string indent)
{
    std::string result = indent + "Node: " + opname + "\n";
    if (!alias.empty())
    {
        result += indent + " (alias: \"" + alias + "\")\n";
    }
    if (!varname.empty())
    {
        result += indent + " (varname: \"" + varname + "\")\n";
    }
    else if (!alias.empty())
    {
        result += indent + " (varname not found)\n";
    }
    if (!indices.empty())
    {
        result += indent + " (indices: [ ";
        for (std::tuple<int, int, int> idx : indices)
        {
            result += std::to_string(std::get<0>(idx)) + ":";
            result += std::to_string(std::get<1>(idx)) + ":";
            result += std::to_string(std::get<2>(idx)) + " ";
        }
        result += "] )\n";
    }
    for (ASTNode *node : sub_exprs)
    {
        result += node->printpretty(indent + "    ");
    }

    return result;
}

std::vector<ASTNode *> ASTNode::get_subexprs() { return sub_exprs; }

std::string ASTNode::get_opname() { return opname; }

std::string ASTNode::get_alias() { return alias; }

std::string ASTNode::get_varname() { return varname; }

std::vector<std::tuple<int, int, int>> ASTNode::get_indices() { return indices; }

double ASTNode::get_value() { return value; }

void ASTNode::set_varname(const std::string s) { varname = s; }

void ASTNode::set_indices(const std::vector<std::tuple<int, int, int>> idx) { indices = idx; }

}
}
