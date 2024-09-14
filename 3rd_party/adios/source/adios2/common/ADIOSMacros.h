/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * ADIOSMacros.h : a set of helper macros used internally
 *
 * Created on: March 23, 2017
 *     Author: Chuck Atkins chuck.atkins@kitware.com
 *             William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_ADIOSMACROS_H
#define ADIOS2_ADIOSMACROS_H

#include <string>

#include "adios2/common/ADIOSTypes.h"

/**
 <pre>
 The ADIOS_FOREACH_TYPE_1ARG macro assumes the given argument is a macro which
 takes a single argument that is a type and then inserts the given MACRO for
 each of the known primitive types

 An example of this might be to instantiate a template function for each
 known type.  For example:

   template<typename T> int foo() { // some implementation of foo  }

   #define instantiate_foo(T) template int foo<T>();
   ADIOS_FOREACH_TYPE_1ARG(instantiate_foo)
   #undef instantiate_foo
 </pre>
*/
#define ADIOS2_FOREACH_ATTRIBUTE_PRIMITIVE_STDTYPE_1ARG(MACRO)                                     \
    MACRO(char)                                                                                    \
    MACRO(int8_t)                                                                                  \
    MACRO(int16_t)                                                                                 \
    MACRO(int32_t)                                                                                 \
    MACRO(int64_t)                                                                                 \
    MACRO(uint8_t)                                                                                 \
    MACRO(uint16_t)                                                                                \
    MACRO(uint32_t)                                                                                \
    MACRO(uint64_t)                                                                                \
    MACRO(float)                                                                                   \
    MACRO(double)                                                                                  \
    MACRO(long double)

#define ADIOS2_FOREACH_PRIMITIVE_STDTYPE_1ARG(MACRO)                                               \
    ADIOS2_FOREACH_ATTRIBUTE_PRIMITIVE_STDTYPE_1ARG(MACRO)                                         \
    MACRO(std::complex<float>)                                                                     \
    MACRO(std::complex<double>)

#define ADIOS2_FOREACH_ATTRIBUTE_STDTYPE_1ARG(MACRO)                                               \
    MACRO(std::string)                                                                             \
    ADIOS2_FOREACH_PRIMITIVE_STDTYPE_1ARG(MACRO)

#define ADIOS2_FOREACH_STDTYPE_1ARG(MACRO)                                                         \
    MACRO(std::string)                                                                             \
    ADIOS2_FOREACH_PRIMITIVE_STDTYPE_1ARG(MACRO)

#define ADIOS2_FOREACH_TYPE_1ARG(MACRO)                                                            \
    MACRO(std::string)                                                                             \
    MACRO(char)                                                                                    \
    MACRO(signed char)                                                                             \
    MACRO(unsigned char)                                                                           \
    MACRO(short)                                                                                   \
    MACRO(unsigned short)                                                                          \
    MACRO(int)                                                                                     \
    MACRO(unsigned int)                                                                            \
    MACRO(long int)                                                                                \
    MACRO(long long int)                                                                           \
    MACRO(unsigned long int)                                                                       \
    MACRO(unsigned long long int)                                                                  \
    MACRO(float)                                                                                   \
    MACRO(double)                                                                                  \
    MACRO(long double)                                                                             \
    MACRO(std::complex<float>)                                                                     \
    MACRO(std::complex<double>)

#define ADIOS2_FOREACH_PRIMITIVE_TYPE_1ARG(MACRO)                                                  \
    MACRO(char)                                                                                    \
    MACRO(signed char)                                                                             \
    MACRO(unsigned char)                                                                           \
    MACRO(short)                                                                                   \
    MACRO(unsigned short)                                                                          \
    MACRO(int)                                                                                     \
    MACRO(unsigned int)                                                                            \
    MACRO(long int)                                                                                \
    MACRO(long long int)                                                                           \
    MACRO(unsigned long int)                                                                       \
    MACRO(unsigned long long int)                                                                  \
    MACRO(float)                                                                                   \
    MACRO(double)                                                                                  \
    MACRO(long double)                                                                             \
    MACRO(std::complex<float>)                                                                     \
    MACRO(std::complex<double>)

#define ADIOS2_FOREACH_COMPLEX_PRIMITIVE_TYPE_1ARG(MACRO)                                          \
    MACRO(float)                                                                                   \
    MACRO(double)                                                                                  \
    MACRO(long double)

#define ADIOS2_FOREACH_CHAR_TYPE_1ARG(MACRO)                                                       \
    MACRO(char)                                                                                    \
    MACRO(signed char)                                                                             \
    MACRO(unsigned char)

#define ADIOS2_FOREACH_NUMERIC_TYPE_1ARG(MACRO)                                                    \
    MACRO(short)                                                                                   \
    MACRO(unsigned short)                                                                          \
    MACRO(int)                                                                                     \
    MACRO(unsigned int)                                                                            \
    MACRO(long int)                                                                                \
    MACRO(long long int)                                                                           \
    MACRO(unsigned long int)                                                                       \
    MACRO(unsigned long long int)                                                                  \
    MACRO(float)                                                                                   \
    MACRO(double)                                                                                  \
    MACRO(long double)                                                                             \
    MACRO(std::complex<float>)                                                                     \
    MACRO(std::complex<double>)

#define ADIOS2_FOREACH_ATTRIBUTE_TYPE_1ARG(MACRO)                                                  \
    MACRO(std::string)                                                                             \
    MACRO(char)                                                                                    \
    MACRO(signed char)                                                                             \
    MACRO(unsigned char)                                                                           \
    MACRO(short)                                                                                   \
    MACRO(unsigned short)                                                                          \
    MACRO(int)                                                                                     \
    MACRO(unsigned int)                                                                            \
    MACRO(long int)                                                                                \
    MACRO(unsigned long int)                                                                       \
    MACRO(long long int)                                                                           \
    MACRO(unsigned long long int)                                                                  \
    MACRO(float)                                                                                   \
    MACRO(double)                                                                                  \
    MACRO(long double)                                                                             \
    MACRO(std::complex<float>)                                                                     \
    MACRO(std::complex<double>)

#define ADIOS2_FOREACH_ATTRIBUTE_PRIMITIVE_TYPE_1ARG(MACRO)                                        \
    MACRO(char)                                                                                    \
    MACRO(signed char)                                                                             \
    MACRO(unsigned char)                                                                           \
    MACRO(short)                                                                                   \
    MACRO(unsigned short)                                                                          \
    MACRO(int)                                                                                     \
    MACRO(unsigned int)                                                                            \
    MACRO(long int)                                                                                \
    MACRO(unsigned long int)                                                                       \
    MACRO(long long int)                                                                           \
    MACRO(unsigned long long int)                                                                  \
    MACRO(float)                                                                                   \
    MACRO(double)                                                                                  \
    MACRO(long double)                                                                             \
    MACRO(std::complex<float>)                                                                     \
    MACRO(std::complex<double>)

/**
 <pre>
 The ADIOS2_FOREACH_STDTYPE_2ARGS macro assumes the given argument is a macro
 which takes two arguments, the first being a type, and the second being a
 label for that type, i.e. std::complex<float> and CFloat.

 An example of this might be to define a virtual method for each known
 type, and since the method is virtual then it cannot be a template.
 For example:

   #define declare_foo(T,L) virtual const T& foo_## L (std::string bar);
   ADIOS2_FOREACH_STDTYPE_2ARGS(declare_foo)
   #undef declare_foo

   is equivalent to:

   virtual           char& foo_Char(std::string bar);
   virtual unsigned  char& foo_UChar(std::string bar);
   virtual          short& foo_Short(std::string bar);
   virtual unsigned short& foo_UShort(std::string bar);
   ...
   virtual std::complex<double>& foo_CDouble(std::string bar);
  </pre>
*/
#define ADIOS2_FOREACH_ATTRIBUTE_STDTYPE_2ARGS(MACRO)                                              \
    MACRO(std::string, string)                                                                     \
    MACRO(int8_t, int8)                                                                            \
    MACRO(char, char)                                                                              \
    MACRO(uint8_t, uint8)                                                                          \
    MACRO(int16_t, int16)                                                                          \
    MACRO(uint16_t, uint16)                                                                        \
    MACRO(int32_t, int32)                                                                          \
    MACRO(uint32_t, uint32)                                                                        \
    MACRO(int64_t, int64)                                                                          \
    MACRO(uint64_t, uint64)                                                                        \
    MACRO(float, float)                                                                            \
    MACRO(double, double)                                                                          \
    MACRO(long double, ldouble)                                                                    \
    MACRO(std::complex<float>, cfloat)                                                             \
    MACRO(std::complex<double>, cdouble)

#define ADIOS2_FOREACH_PRIMITVE_STDTYPE_2ARGS(MACRO)                                               \
    MACRO(int8_t, int8)                                                                            \
    MACRO(char, char)                                                                              \
    MACRO(uint8_t, uint8)                                                                          \
    MACRO(int16_t, int16)                                                                          \
    MACRO(uint16_t, uint16)                                                                        \
    MACRO(int32_t, int32)                                                                          \
    MACRO(uint32_t, uint32)                                                                        \
    MACRO(int64_t, int64)                                                                          \
    MACRO(uint64_t, uint64)                                                                        \
    MACRO(float, float)                                                                            \
    MACRO(double, double)                                                                          \
    MACRO(long double, ldouble)                                                                    \
    MACRO(std::complex<float>, cfloat)                                                             \
    MACRO(std::complex<double>, cdouble)

#define ADIOS2_FOREACH_MINMAX_STDTYPE_2ARGS(MACRO)                                                 \
    MACRO(int8_t, int8)                                                                            \
    MACRO(uint8_t, uint8)                                                                          \
    MACRO(int16_t, int16)                                                                          \
    MACRO(uint16_t, uint16)                                                                        \
    MACRO(int32_t, int32)                                                                          \
    MACRO(uint32_t, uint32)                                                                        \
    MACRO(int64_t, int64)                                                                          \
    MACRO(uint64_t, uint64)                                                                        \
    MACRO(float, float)                                                                            \
    MACRO(double, double)                                                                          \
    MACRO(long double, ldouble)

#define ADIOS2_FOREACH_STDTYPE_2ARGS(MACRO) ADIOS2_FOREACH_ATTRIBUTE_STDTYPE_2ARGS(MACRO)

#define ADIOS2_CLASS_iterator                                                                      \
    class iterator                                                                                 \
    {                                                                                              \
    public:                                                                                        \
        typedef iterator self_type;                                                                \
        typedef T value_type;                                                                      \
        typedef T &reference;                                                                      \
        typedef T *pointer;                                                                        \
        typedef std::forward_iterator_tag iterator_category;                                       \
        typedef int difference_type;                                                               \
        iterator(pointer ptr) : ptr_(ptr) {}                                                       \
        self_type operator++()                                                                     \
        {                                                                                          \
            self_type i = *this;                                                                   \
            ptr_++;                                                                                \
            return i;                                                                              \
        }                                                                                          \
        self_type operator++(int junk)                                                             \
        {                                                                                          \
            ptr_++;                                                                                \
            return *this;                                                                          \
        }                                                                                          \
        reference operator*() { return *ptr_; }                                                    \
        pointer operator->() { return ptr_; }                                                      \
        bool operator==(const self_type &rhs) { return ptr_ == rhs.ptr_; }                         \
        bool operator!=(const self_type &rhs) { return ptr_ != rhs.ptr_; }                         \
                                                                                                   \
    private:                                                                                       \
        pointer ptr_;                                                                              \
    }

#define ADIOS2_iterators_functions(DATA_FUNCTION, SIZE_FUNCTION)                                   \
    iterator begin() noexcept { return iterator(DATA_FUNCTION); }                                  \
    iterator end() noexcept { return iterator(DATA_FUNCTION + SIZE_FUNCTION); }

#endif /* ADIOS2_ADIOSMACROS_H */
