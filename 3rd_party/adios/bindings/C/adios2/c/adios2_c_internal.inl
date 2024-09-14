
/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adios2_c_internal.inl
 *
 *  Created on: Feb 9, 2019
 *      Author: Kai Germaschewski <kai.germaschewski@unh.edu>
 */

#ifndef ADIOS2_BINDINGS_C_C_ADIOS2_C_INTERNAL_INL_
#define ADIOS2_BINDINGS_C_C_ADIOS2_C_INTERNAL_INL_
#ifndef ADIOS2_BINDINGS_C_C_ADIOS2_C_INTERNAL_H_
#error "Inline file should only be included from its header, never on its own"
#endif

#include <complex>

namespace
{

// no default implementation, so only type below are
// supported

#define make_MapAdios2Type(adios2_type, T)                                     \
    template <>                                                                \
    struct MapAdios2Type<adios2_type>                                          \
    {                                                                          \
        using Type = T;                                                        \
    };
/* clang-format off */
make_MapAdios2Type(adios2_type_int8_t, int8_t)
make_MapAdios2Type(adios2_type_int16_t, int16_t)
make_MapAdios2Type(adios2_type_int32_t, int32_t)
make_MapAdios2Type(adios2_type_int64_t, int64_t)
make_MapAdios2Type(adios2_type_uint8_t, uint8_t)
make_MapAdios2Type(adios2_type_uint16_t, uint16_t)
make_MapAdios2Type(adios2_type_uint32_t, uint32_t)
make_MapAdios2Type(adios2_type_uint64_t, uint64_t)
make_MapAdios2Type(adios2_type_float, float)
make_MapAdios2Type(adios2_type_double, double)
make_MapAdios2Type(adios2_type_long_double, long double)
make_MapAdios2Type(adios2_type_float_complex, std::complex<float>)
make_MapAdios2Type(adios2_type_double_complex, std::complex<double>)
/* clang-format on */
#undef make_MapAdios2Type

                                    inline adios2_error
    String2CAPI(const std::string &s, char *buf, size_t *size)
{
    *size = s.size();
    if (buf != nullptr)
    {
        s.copy(buf, *size);
    }
    return adios2_error_none;
}

} // namespace

#endif /* ADIOS2_BINDINGS_C_C_ADIOS2_C_INTERNAL_INL_ */
