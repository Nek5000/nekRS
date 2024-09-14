#include "ASTDriver.h"

namespace adios2
{
namespace detail
{

using indx_type = std::vector<std::tuple<int, int, int>>;

ASTDriver::ASTDriver() {}

ASTDriver::ASTDriver(const std::string input) { ASTDriver::parse(input); }

ASTDriver::~ASTDriver()
{
    while (holding.size() > 0)
    {
        delete holding.top();
        holding.pop();
    }
}

ASTNode *ASTDriver::getAST()
{
    // TODO: check only one ASTNode remains in holding
    // else throw error that parsing failed
    resolve(holding.top());
    return holding.top();
}

void ASTDriver::resolve(ASTNode *node)
{
    if (!node->get_alias().empty())
    {
        std::tuple<std::string, indx_type> var_info;
        var_info = lookup_var(node->get_alias());
        node->set_varname(std::get<0>(var_info));
        node->set_indices(std::get<1>(var_info));
    }
    for (ASTNode *subexpr : node->get_subexprs())
    {
        resolve(subexpr);
    }
}

std::tuple<std::string, indx_type> ASTDriver::lookup_var(const std::string alias)
{
    return aliases[alias];
}

std::string ASTDriver::lookup_var_name(const std::string alias)
{
    std::tuple<std::string, indx_type> var = aliases[alias];
    return std::get<0>(var);
}

indx_type ASTDriver::lookup_var_indices(const std::string alias)
{
    std::tuple<std::string, indx_type> var = aliases[alias];
    return std::get<1>(var);
}

void ASTDriver::add_lookup_entry(std::string alias, std::string var_name, indx_type indices)
{
    aliases.insert({alias, {var_name, indices}});
}

void ASTDriver::add_lookup_entry(std::string alias, std::string var_name)
{
    aliases.insert({alias, {var_name, {}}});
}

void ASTDriver::createNode(std::string op_name, size_t numsubexprs)
{
    ASTNode *node = new ASTNode(op_name, numsubexprs);
    for (size_t i = 1; i <= numsubexprs; ++i)
    {
        // TODO: check that holding contains ASTNode(s)
        // else throw error that parsing failed
        ASTNode *subexpr = holding.top();
        node->insert_subexpr_n(subexpr, numsubexprs - i);
        holding.pop();
    }
    holding.push(node);
}

void ASTDriver::createNode(std::string alias)
{
    ASTNode *node = new ASTNode("ALIAS", alias);
    holding.push(node);
}

void ASTDriver::createNode(std::string alias, indx_type indices)
{
    ASTNode *node = new ASTNode("INDEX", indices);
    node->pushback_subexpr(new ASTNode("ALIAS", alias));
    holding.push(node);
}

}
}
