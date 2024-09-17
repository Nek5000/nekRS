#ifndef ADIOS2_DERIVED_PARSER_EXPHELPER_H_
#define ADIOS2_DERIVED_PARSER_EXPHELPER_H_

#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace adios2
{
namespace detail
{

enum ExpressionOperator
{
    OP_NULL,
    OP_ALIAS, /* Parser-use only */
    OP_PATH,  /* Parser-use only */
    OP_NUM,   /* Parser-use only */
    OP_INDEX,
    OP_ADD,
    OP_SQRT,
    OP_POW,
    OP_MAGN,
    OP_CURL
};

struct OperatorProperty
{
    std::string name;
    bool is_associative;
};

const std::map<ExpressionOperator, OperatorProperty> op_property = {
    {ExpressionOperator::OP_NULL, {"NULL", false}},
    {ExpressionOperator::OP_ALIAS, {"ALIAS", false}}, /* Parser-use only */
    {ExpressionOperator::OP_PATH, {"PATH", false}},   /* Parser-use only */
    {ExpressionOperator::OP_NUM, {"NUM", false}},     /* Parser-use only */
    {ExpressionOperator::OP_INDEX, {"INDEX", false}},
    {ExpressionOperator::OP_ADD, {"ADD", true}},
    {ExpressionOperator::OP_SQRT, {"SQRT", false}},
    {ExpressionOperator::OP_POW, {"POW", false}},
    {ExpressionOperator::OP_CURL, {"CURL", false}},
    {ExpressionOperator::OP_MAGN, {"MAGNITUDE", false}}};

const std::map<std::string, ExpressionOperator> string_to_op = {
    {"ALIAS", ExpressionOperator::OP_ALIAS}, /* Parser-use only */
    {"PATH", ExpressionOperator::OP_PATH},   /* Parser-use only */
    {"NUM", ExpressionOperator::OP_NUM},     /* Parser-use only */
    {"INDEX", ExpressionOperator::OP_INDEX},    {"+", ExpressionOperator::OP_ADD},
    {"add", ExpressionOperator::OP_ADD},        {"ADD", ExpressionOperator::OP_ADD},
    {"SQRT", ExpressionOperator::OP_SQRT},      {"sqrt", ExpressionOperator::OP_SQRT},
    {"POW", ExpressionOperator::OP_POW},        {"^", ExpressionOperator::OP_POW},
    {"CURL", ExpressionOperator::OP_CURL},      {"curl", ExpressionOperator::OP_CURL},
    {"MAGNITUDE", ExpressionOperator::OP_MAGN}, {"magnitude", ExpressionOperator::OP_MAGN}};

inline std::string get_op_name(ExpressionOperator op) { return op_property.at(op).name; }

inline ExpressionOperator get_op(std::string op) { return string_to_op.at(op); }

}
}
#endif
