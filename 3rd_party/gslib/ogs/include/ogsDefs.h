/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

/* the supported types */
typedef long long long_long;
#if 0
#define OGS_FOR_EACH_TYPE(macro) \
  macro(double   ) \
  macro(float    ) \
  macro(int      ) \
  macro(long     ) \
  macro(long_long)
#else
#define OGS_FOR_EACH_TYPE(macro) \
  macro(double   ) \
  macro(float    )
#endif

/* the supported ops */
#if 0
#define OGS_FOR_EACH_OP(T,macro) \
  macro(T,add) \
  macro(T,mul) \
  macro(T,min) \
  macro(T,max)
#else
#define OGS_FOR_EACH_OP(T,macro) \
  macro(T,add) \
  macro(T,min) \
  macro(T,max)
#endif

#define OGS_DO_add(a,b) a+=b
#define OGS_DO_mul(a,b) a*=b
#define OGS_DO_min(a,b) if(b<a) a=b
#define OGS_DO_max(a,b) if(b>a) a=b

/* type size array */
#define OGS_TYPE_SIZE_ITEM(T) sizeof(T),
#define OGS_DEFINE_TYPE_SIZES() \
  static const unsigned ogs_type_size[] = \
    { OGS_FOR_EACH_TYPE(OGS_TYPE_SIZE_ITEM) 0 };

/* mapping from ogs types to gs types */
#define gs_int64_t gs_long_long
#define OGS_GS_MAP_TYPE_ITEM(T) gs_##T,
#define OGS_GS_DEFINE_TYPE_MAP() \
  static const gs_dom ogs_gs_type_map[] = \
    { OGS_FOR_EACH_TYPE(OGS_GS_MAP_TYPE_ITEM) gs_dom_n };

/* mapping from ogs ops to gs ops */
#define OGS_GS_MAP_OP_ITEM(T,OP) gs_##OP,
#define OGS_GS_DEFINE_OP_MAP() \
  static const gs_op ogs_gs_op_map[] = \
    { OGS_FOR_EACH_OP(T,OGS_GS_MAP_OP_ITEM) gs_op_n };
