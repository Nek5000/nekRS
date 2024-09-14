#ifndef ASTDRIVER_HH_
#define ASTDRIVER_HH_

#include "ASTNode.h"
#include "parser.h"
#include <map>
#include <stack>
#include <string>
#include <tuple>

#define YY_DECL adios2::detail::parser::symbol_type yylex(adios2::detail::ASTDriver &drv)
YY_DECL;

namespace adios2
{
namespace detail
{

using indx_type = std::vector<std::tuple<int, int, int>>;

class ASTDriver
{
public:
    ASTDriver();
    ASTDriver(const std::string input);
    ~ASTDriver();

    // Defined in lexer.l
    void parse(const std::string input);

    ASTNode *getAST();

    void resolve(ASTNode *node);

    std::tuple<std::string, indx_type> lookup_var(const std::string alias);
    std::string lookup_var_name(const std::string alias);
    indx_type lookup_var_indices(const std::string alias);

    void add_lookup_entry(std::string alias, std::string var_name, indx_type indices);
    void add_lookup_entry(std::string alias, std::string var_name);

    void createNode(std::string, size_t);
    void createNode(std::string);
    void createNode(std::string, indx_type);

    // Whether to generate parser debug traces.
    bool trace_parsing = false;
    // Whether to generate scanner debug traces.
    bool trace_scanning = false;
    // The token's location used by the scanner.
    adios2::detail::location location;

private:
    // While parsing, holds ASTNodes until parent node is created
    // (since root node is created last from bottom up design)
    std::stack<ASTNode *> holding;

    // Variable lookup table: maps alias names
    // to variable names and indices from alias definition
    std::map<std::string, std::tuple<std::string, indx_type>> aliases;
};

}
}
#endif // ! ASTDRIVER_HH_
