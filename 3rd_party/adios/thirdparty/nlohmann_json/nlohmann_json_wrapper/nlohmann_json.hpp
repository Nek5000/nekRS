#ifndef __NLOHMANN_JSON_HPP__
#define __NLOHMANN_JSON_HPP__

#if defined(_MSC_VER) && !defined(__clang__)
#pragma warning(push)

#elif defined(__INTEL_COMPILER) && !defined(__clang__)
#pragma warning push
#pragma warning disable 1011
_Pragma("warning(disable:1011)")

#elif defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"

#elif defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wall"
#endif

#include <nlohmann/json.hpp>

#if defined(_MSC_VER) && !defined(__clang__)
#pragma warning(pop)

#elif defined(__INTEL_COMPILER) && !defined(__clang__)
#pragma warning pop

#elif defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop

#elif defined(__clang__)
#pragma clang diagnostic pop
#endif

#endif
