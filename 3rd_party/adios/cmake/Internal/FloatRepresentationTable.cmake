# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

# The intention of these magic values is to find the boundaries between the
# different components of a floating point number. There are four finite
# numbers: one baseline, one sign inversion, one exponent inversion, and
# one significand change (not a perfect inversion). The sign is positive
# (0) in the baseline, and negative (1) in the inversion. The exponent is 0
# (which evaluates to 0[1...]) in the baseline and 1 (which evaluates to
# 1[0...]) in the inversion. The significand is 1.025 (which evaluates to
# 10000[0110...]) in the baseline and 1.9875 (which evaluates to
# 11111[1100...]) in the inversion. This value was chosen for the
# significand because it ensures that the least significant bit is always
# inverted no matter where it gets rounded off.
#
# Testing NaN is tricky, because the significand bits can be any non-zero
# value (and compilers do emit different values in practice). While this can
# theoretically be tested with a regular expression, in practice such a regex
# would be exceedingly long, so we simply test for any valid significand and
# hope that it's different from INFINITY.
#
# Note: We will skip the INFINITY and NAN tests for now because there is not a
# portable way to test these on Fortran.
set(_tests
  BASELINE
  SIGN
  EXPONENT
  SIGNIFICAND
  #INFINITY
  #N_INFINITY
  #NAN
  )
if(_cftr_LANGUAGE STREQUAL "Fortran")
  set(_test_BASELINE
    "1.025q0"
    "1.025dd0"
    "1.025d0"
    "1.025s0"
    )
  set(_test_SIGN
    "-1.025q0"
    "-1.025dd0"
    "-1.025d0"
    "-1.025s0"
    )
  set(_test_EXPONENT
    "2.05q0"
    "2.05dd0"
    "2.05d0"
    "2.05s0"
    )
  set(_test_SIGNIFICAND
    "1.9875q0"
    "1.9875dd0"
    "1.9875d0"
    "1.9875s0"
    )
  #set(_test_INFINITY    "INFINITY")
  #set(_test_N_INFINITY  "-INFINITY")
  #set(_test_NAN         "NAN")
else()
  set(_test_BASELINE    "1.025l")
  set(_test_SIGN        "-1.025l")
  set(_test_EXPONENT    "2.05l")
  set(_test_SIGNIFICAND "1.9875l")
  #set(_test_INFINITY    "INFINITY")
  #set(_test_N_INFINITY  "-INFINITY")
  #set(_test_NAN         "NAN")
endif()

set(_reprs
  FLOAT_HP_B_IEEE754_LE
  FLOAT_HP_B_IEEE754_BE
  FLOAT_SP_B_IEEE754_LE
  FLOAT_SP_B_IEEE754_BE
  FLOAT_DP_B_IEEE754_LE
  FLOAT_DP_B_IEEE754_BE
  FLOAT_QP_B_IEEE754_LE
  FLOAT_QP_B_IEEE754_BE
  FLOAT_QP_B_IEEE754_80_LE
  FLOAT_EP_B_X86_LE
  FLOAT_EP_B_X86_64_LE
  FLOAT_EP_B_IBM_EXTENDED_LE
  FLOAT_EP_B_IBM_EXTENDED_BE
  )

# IEEE 754 half precision - 5 exponent bits, 11 significand bits (10 explicit)
set(_repr_FLOAT_HP_B_IEEE754_LE_BASELINE    "^1a3c$")
set(_repr_FLOAT_HP_B_IEEE754_LE_SIGN        "^1abc$")
set(_repr_FLOAT_HP_B_IEEE754_LE_EXPONENT    "^1a40$")
set(_repr_FLOAT_HP_B_IEEE754_LE_SIGNIFICAND "^f33f$")
#set(_repr_FLOAT_HP_B_IEEE754_LE_INFINITY    "^007c$")
#set(_repr_FLOAT_HP_B_IEEE754_LE_N_INFINITY  "^00fc$")
#set(_repr_FLOAT_HP_B_IEEE754_LE_NAN         "^[0-9a-f]+7[c-f]$")

set(_repr_FLOAT_HP_B_IEEE754_BE_BASELINE    "^3c1a$")
set(_repr_FLOAT_HP_B_IEEE754_BE_SIGN        "^bc1a$")
set(_repr_FLOAT_HP_B_IEEE754_BE_EXPONENT    "^401a$")
set(_repr_FLOAT_HP_B_IEEE754_BE_SIGNIFICAND "^3ff3$")
#set(_repr_FLOAT_HP_B_IEEE754_BE_INFINITY    "^7c00$")
#set(_repr_FLOAT_HP_B_IEEE754_BE_N_INFINITY  "^fc00$")
#set(_repr_FLOAT_HP_B_IEEE754_BE_NAN         "^7[c-f][0-9a-f]+$")

# IEEE 754 single precision - 8 exponent bits, 24 significand bits (23 explicit)
set(_repr_FLOAT_SP_B_IEEE754_LE_BASELINE    "^3333833f$")
set(_repr_FLOAT_SP_B_IEEE754_LE_SIGN        "^333383bf$")
set(_repr_FLOAT_SP_B_IEEE754_LE_EXPONENT    "^33330340$")
set(_repr_FLOAT_SP_B_IEEE754_LE_SIGNIFICAND "^6666fe3f$")
#set(_repr_FLOAT_SP_B_IEEE754_LE_INFINITY    "^0000807f$")
#set(_repr_FLOAT_SP_B_IEEE754_LE_N_INFINITY  "^000080ff$")
#set(_repr_FLOAT_SP_B_IEEE754_LE_NAN         "^[0-9a-f]+[89a-f][0-9a-f]7f$")

set(_repr_FLOAT_SP_B_IEEE754_BE_BASELINE    "^3f833333$")
set(_repr_FLOAT_SP_B_IEEE754_BE_SIGN        "^bf833333$")
set(_repr_FLOAT_SP_B_IEEE754_BE_EXPONENT    "^40033333$")
set(_repr_FLOAT_SP_B_IEEE754_BE_SIGNIFICAND "^3ffe6666$")
#set(_repr_FLOAT_SP_B_IEEE754_BE_INFINITY    "^7f800000$")
#set(_repr_FLOAT_SP_B_IEEE754_BE_N_INFINITY  "^ff800000$")
#set(_repr_FLOAT_SP_B_IEEE754_BE_NAN         "^7f[89a-f][0-9a-f]+$")

# IEEE 754 double precision - 11 exponent bits, 53 significand bits (52 explicit)
set(_repr_FLOAT_DP_B_IEEE754_LE_BASELINE    "^666666666666f03f$")
set(_repr_FLOAT_DP_B_IEEE754_LE_SIGN        "^666666666666f0bf$")
set(_repr_FLOAT_DP_B_IEEE754_LE_EXPONENT    "^6666666666660040$")
set(_repr_FLOAT_DP_B_IEEE754_LE_SIGNIFICAND "^cdccccccccccff3f$")
#set(_repr_FLOAT_DP_B_IEEE754_LE_INFINITY    "^000000000000f07f$")
#set(_repr_FLOAT_DP_B_IEEE754_LE_N_INFINITY  "^000000000000f0ff$")
#set(_repr_FLOAT_DP_B_IEEE754_LE_NAN         "^[0-9a-f]+f[0-9a-f]7f$")

set(_repr_FLOAT_DP_B_IEEE754_BE_BASELINE    "^3ff0666666666666$")
set(_repr_FLOAT_DP_B_IEEE754_BE_SIGN        "^bff0666666666666$")
set(_repr_FLOAT_DP_B_IEEE754_BE_EXPONENT    "^4000666666666666$")
set(_repr_FLOAT_DP_B_IEEE754_BE_SIGNIFICAND "^3fffcccccccccccd$")
#set(_repr_FLOAT_DP_B_IEEE754_BE_INFINITY    "^7ff0000000000000$")
#set(_repr_FLOAT_DP_B_IEEE754_BE_N_INFINITY  "^fff0000000000000$")
#set(_repr_FLOAT_DP_B_IEEE754_BE_NAN         "^7ff[0-9a-f]+$")

# IEEE 754 quadruple precision - 15 exponent bits, 113 significand bits (112 explicit)
set(_repr_FLOAT_QP_B_IEEE754_LE_BASELINE    "^6666666666666666666666666606ff3f$")
set(_repr_FLOAT_QP_B_IEEE754_LE_SIGN        "^6666666666666666666666666606ffbf$")
set(_repr_FLOAT_QP_B_IEEE754_LE_EXPONENT    "^66666666666666666666666666060040$")
set(_repr_FLOAT_QP_B_IEEE754_LE_SIGNIFICAND "^cdccccccccccccccccccccccccfcff3f$")
#set(_repr_FLOAT_QP_B_IEEE754_LE_INFINITY    "^0000000000000000000000000000ff7f$")
#set(_repr_FLOAT_QP_B_IEEE754_LE_N_INFINITY  "^0000000000000000000000000000ffff$")
#set(_repr_FLOAT_QP_B_IEEE754_LE_NAN         "^[0-9a-f]+ff7f$")

set(_repr_FLOAT_QP_B_IEEE754_BE_BASELINE    "^3fff0666666666666666666666666666$")
set(_repr_FLOAT_QP_B_IEEE754_BE_SIGN        "^bfff0666666666666666666666666666$")
set(_repr_FLOAT_QP_B_IEEE754_BE_EXPONENT    "^40000666666666666666666666666666$")
set(_repr_FLOAT_QP_B_IEEE754_BE_SIGNIFICAND "^3ffffccccccccccccccccccccccccccd$")
#set(_repr_FLOAT_QP_B_IEEE754_BE_INFINITY    "^7fff0000000000000000000000000000$")
#set(_repr_FLOAT_QP_B_IEEE754_BE_N_INFINITY  "^ffff0000000000000000000000000000$")
#set(_repr_FLOAT_QP_B_IEEE754_BE_NAN         "^7fff[0-9a-f]+$")

# IEEE 754 quadruple precision with only 80 bits
set(_repr_FLOAT_QP_B_IEEE754_80_LE_BASELINE     "^0000000000006666666666666606ff3f$")
set(_repr_FLOAT_QP_B_IEEE754_80_LE_SIGN         "^0000000000006666666666666606ffbf$")
set(_repr_FLOAT_QP_B_IEEE754_80_LE_EXPONENT     "^00000000000066666666666666060040$")
set(_repr_FLOAT_QP_B_IEEE754_80_LE_SIGNIFICAND  "^000000000000ccccccccccccccfcff3f$")
#set(_repr_FLOAT_QP_B_IEEE754_80_LE_INFINITY     "^0000000000000000000000000000ff7f$")
#set(_repr_FLOAT_QP_B_IEEE754_80_LE_N_INFINITY   "^0000000000000000000000000000ffff$")
#set(_repr_FLOAT_QP_B_IEEE754_80_LE_NAN          "^[0-9a-f]+ff7f$")

# x86 extended precision - 15 exponent bits, 64 significand bits (all explicit)
set(_repr_FLOAT_EP_B_X86_LE_BASELINE    "^3333333333333383ff3f0000$")
set(_repr_FLOAT_EP_B_X86_LE_SIGN        "^3333333333333383ffbf0000$")
set(_repr_FLOAT_EP_B_X86_LE_EXPONENT    "^333333333333338300400000$")
set(_repr_FLOAT_EP_B_X86_LE_SIGNIFICAND "^66666666666666feff3f0000$")
#set(_repr_FLOAT_EP_B_X86_LE_INFINITY    "^0000000000000080ff7f0000$")
#set(_repr_FLOAT_EP_B_X86_LE_N_INFINITY  "^0000000000000080ffff0000$")
#set(_repr_FLOAT_EP_B_X86_LE_NAN         "^[0-9a-f]+[89a-f][0-9a-f]ff7f0000$")

set(_repr_FLOAT_EP_B_X86_64_LE_BASELINE     "^3333333333333383ff3f000000000000$")
set(_repr_FLOAT_EP_B_X86_64_LE_SIGN         "^3333333333333383ffbf000000000000$")
set(_repr_FLOAT_EP_B_X86_64_LE_EXPONENT     "^33333333333333830040000000000000$")
set(_repr_FLOAT_EP_B_X86_64_LE_SIGNIFICAND  "^66666666666666feff3f000000000000$")
#set(_repr_FLOAT_EP_B_X86_64_LE_INFINITY     "^0000000000000080ff7f000000000000$")
#set(_repr_FLOAT_EP_B_X86_64_LE_N_INFINITY   "^0000000000000080ffff000000000000$")
#set(_repr_FLOAT_EP_B_X86_64_LE_NAN          "^[0-9a-f]+[89a-f][0-9a-f]ff7f000000000000$")

# IBM extended precision - two IEEE 754 doubles back to back
set(_repr_FLOAT_EP_B_IBM_EXTENDED_LE_BASELINE     "^666666666666f03f9a9999999999993c$")
set(_repr_FLOAT_EP_B_IBM_EXTENDED_LE_SIGN         "^666666666666f0bf9a999999999999bc$")
set(_repr_FLOAT_EP_B_IBM_EXTENDED_LE_EXPONENT     "^66666666666600409a9999999999a93c$")
set(_repr_FLOAT_EP_B_IBM_EXTENDED_LE_SIGNIFICAND  "^cdccccccccccff3f9[8a]999999999989bc$")
#set(_repr_FLOAT_EP_B_IBM_EXTENDED_LE_INFINITY     "^000000000000f07f0000000000000000$")
#set(_repr_FLOAT_EP_B_IBM_EXTENDED_LE_N_INFINITY   "^000000000000f0ff0000000000000000$")
#set(_repr_FLOAT_EP_B_IBM_EXTENDED_LE_NAN          "^[0-9a-f]+f[0-9a-f]7f[0-9a-f]+$")

set(_repr_FLOAT_EP_B_IBM_EXTENDED_BE_BASELINE     "^3ff06666666666663c9999999999999a$")
set(_repr_FLOAT_EP_B_IBM_EXTENDED_BE_SIGN         "^bff0666666666666bc9999999999999a$")
set(_repr_FLOAT_EP_B_IBM_EXTENDED_BE_EXPONENT     "^40006666666666663ca999999999999a$")
set(_repr_FLOAT_EP_B_IBM_EXTENDED_BE_SIGNIFICAND  "^3fffcccccccccccdbc89999999999998$")
#set(_repr_FLOAT_EP_B_IBM_EXTENDED_BE_INFINITY     "^7ff00000000000000000000000000000$")
#set(_repr_FLOAT_EP_B_IBM_EXTENDED_BE_N_INFINITY   "^fff00000000000000000000000000000$")
#set(_repr_FLOAT_EP_B_IBM_EXTENDED_BE_NAN          "^7ff[0-9a-f]+$")
