
/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * ADIOSTypes.inl : inline implementatios for ADIOSTypes.h
 *
 *  Created on: Feb 20, 2019
 *      Author: Kai Germaschewski <kai.germaschewski@unh.edu>
 *              Chuck Atkins chuck.atkins@kitware.com
 *              Norbert Podhorszki pnorbert@ornl.gov
 *              William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_ADIOSTYPES_INL_
#define ADIOS2_ADIOSTYPES_INL_
#ifndef ADIOS2_ADIOSTYPES_H_
#error "Inline file should only be included from its header, never on its own"
#endif

namespace adios2
{

namespace
{

// Get a fixed width integer type from a size specification
template <size_t Bytes, bool Signed>
struct FixedWidthInt;

template <>
struct FixedWidthInt<1, true>
{
    using Type = std::int8_t;
};
template <>
struct FixedWidthInt<2, true>
{
    using Type = std::int16_t;
};
template <>
struct FixedWidthInt<4, true>
{
    using Type = std::int32_t;
};
template <>
struct FixedWidthInt<8, true>
{
    using Type = std::int64_t;
};
template <>
struct FixedWidthInt<1, false>
{
    using Type = std::uint8_t;
};
template <>
struct FixedWidthInt<2, false>
{
    using Type = std::uint16_t;
};
template <>
struct FixedWidthInt<4, false>
{
    using Type = std::uint32_t;
};
template <>
struct FixedWidthInt<8, false>
{
    using Type = std::uint64_t;
};

} // namespace

// Some core type information that may be useful at compile time
template <typename T, typename Enable>
struct TypeInfo
{
    using IOType = T;
    using ValueType = T;
};

// Hack "char" type into this convoluted struct definition which is the key
// to translate between bindings API types to supported stdtypes
template <>
struct TypeInfo<char,
                typename std::enable_if<std::is_same<char, char>::value>::type>
{
    using IOType = char;
    using ValueType = char;
};

template <typename T>
struct TypeInfo<T, typename std::enable_if<std::is_integral<T>::value>::type>
{
    using IOType =
        typename FixedWidthInt<sizeof(T), std::is_signed<T>::value>::Type;
    using ValueType = T;
};

template <typename T>
struct TypeInfo<T,
                typename std::enable_if<std::is_floating_point<T>::value>::type>
{
    using IOType = T;
    using ValueType = T;
};

template <typename T>
struct TypeInfo<T, typename std::enable_if<std::is_same<
                       T, std::complex<typename T::value_type>>::value>::type>
{
    using IOType = T;
    using ValueType = typename T::value_type;
};

template <typename T>
struct TypeInfo<
    T, typename std::enable_if<std::is_same<T, std::string>::value>::type>
{
    using IOType = T;
    using ValueType = T;
};

template <typename T, typename Enable>
inline std::ostream &operator<<(std::ostream &os, const T &value)
{
    return os << ToString(value);
}

} // end namespace adios2

#endif /* ADIOS2_ADIOSTYPES_INL_ */
