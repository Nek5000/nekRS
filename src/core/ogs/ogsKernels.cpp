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

#include "ogstypes.h"
#include "ogs.hpp"
#include "ogsKernels.hpp"
#include "platform.hpp"

namespace ogs {

int Nrefs = 0;
const int gatherNodesPerBlock =
    std::min(BLOCKSIZE, 1024); // needs to be tuned based on arch

void *hostBuf;
size_t hostBufSize = 0;

void *haloBuf;
occa::memory o_haloBuf;
occa::memory h_haloBuf;

occa::stream defaultStream;
occa::stream dataStream;

occa::properties kernelInfo;

occa::kernel gatherScatterNewKernel_floatAdd;
occa::kernel gatherScatterNewKernel_doubleAdd;
occa::kernel gatherScatterNewKernel_doubleMin;
occa::kernel gatherScatterNewKernel_doubleMax;

occa::kernel gatherScatterKernel_floatAdd;
occa::kernel gatherScatterKernel_floatMul;
occa::kernel gatherScatterKernel_floatMin;
occa::kernel gatherScatterKernel_floatMax;
occa::kernel gatherScatterKernel_doubleAdd;
occa::kernel gatherScatterKernel_doubleMul;
occa::kernel gatherScatterKernel_doubleMin;
occa::kernel gatherScatterKernel_doubleMax;
occa::kernel gatherScatterKernel_intAdd;
occa::kernel gatherScatterKernel_intMul;
occa::kernel gatherScatterKernel_intMin;
occa::kernel gatherScatterKernel_intMax;
occa::kernel gatherScatterKernel_longAdd;
occa::kernel gatherScatterKernel_longMul;
occa::kernel gatherScatterKernel_longMin;
occa::kernel gatherScatterKernel_longMax;
occa::kernel gatherScatterVecKernel_floatAdd;
occa::kernel gatherScatterVecKernel_floatMul;
occa::kernel gatherScatterVecKernel_floatMin;
occa::kernel gatherScatterVecKernel_floatMax;
occa::kernel gatherScatterVecKernel_doubleAdd;
occa::kernel gatherScatterVecKernel_doubleMul;
occa::kernel gatherScatterVecKernel_doubleMin;
occa::kernel gatherScatterVecKernel_doubleMax;
occa::kernel gatherScatterVecKernel_intAdd;
occa::kernel gatherScatterVecKernel_intMul;
occa::kernel gatherScatterVecKernel_intMin;
occa::kernel gatherScatterVecKernel_intMax;
occa::kernel gatherScatterVecKernel_longAdd;
occa::kernel gatherScatterVecKernel_longMul;
occa::kernel gatherScatterVecKernel_longMin;
occa::kernel gatherScatterVecKernel_longMax;
occa::kernel gatherScatterManyKernel_floatAdd;
occa::kernel gatherScatterManyKernel_floatMul;
occa::kernel gatherScatterManyKernel_floatMin;
occa::kernel gatherScatterManyKernel_floatMax;
occa::kernel gatherScatterManyKernel_doubleAdd;
occa::kernel gatherScatterManyKernel_doubleMul;
occa::kernel gatherScatterManyKernel_doubleMin;
occa::kernel gatherScatterManyKernel_doubleMax;
occa::kernel gatherScatterManyKernel_intAdd;
occa::kernel gatherScatterManyKernel_intMul;
occa::kernel gatherScatterManyKernel_intMin;
occa::kernel gatherScatterManyKernel_intMax;
occa::kernel gatherScatterManyKernel_longAdd;
occa::kernel gatherScatterManyKernel_longMul;
occa::kernel gatherScatterManyKernel_longMin;
occa::kernel gatherScatterManyKernel_longMax;

occa::kernel gatherKernel_floatAdd;
occa::kernel gatherKernel_floatMul;
occa::kernel gatherKernel_floatMin;
occa::kernel gatherKernel_floatMax;
occa::kernel gatherKernel_doubleAdd;
occa::kernel gatherKernel_doubleMul;
occa::kernel gatherKernel_doubleMin;
occa::kernel gatherKernel_doubleMax;
occa::kernel gatherKernel_intAdd;
occa::kernel gatherKernel_intMul;
occa::kernel gatherKernel_intMin;
occa::kernel gatherKernel_intMax;
occa::kernel gatherKernel_longAdd;
occa::kernel gatherKernel_longMul;
occa::kernel gatherKernel_longMin;
occa::kernel gatherKernel_longMax;
occa::kernel gatherVecKernel_floatAdd;
occa::kernel gatherVecKernel_floatMul;
occa::kernel gatherVecKernel_floatMin;
occa::kernel gatherVecKernel_floatMax;
occa::kernel gatherVecKernel_doubleAdd;
occa::kernel gatherVecKernel_doubleMul;
occa::kernel gatherVecKernel_doubleMin;
occa::kernel gatherVecKernel_doubleMax;
occa::kernel gatherVecKernel_intAdd;
occa::kernel gatherVecKernel_intMul;
occa::kernel gatherVecKernel_intMin;
occa::kernel gatherVecKernel_intMax;
occa::kernel gatherVecKernel_longAdd;
occa::kernel gatherVecKernel_longMul;
occa::kernel gatherVecKernel_longMin;
occa::kernel gatherVecKernel_longMax;
occa::kernel gatherManyKernel_floatAdd;
occa::kernel gatherManyKernel_floatMul;
occa::kernel gatherManyKernel_floatMin;
occa::kernel gatherManyKernel_floatMax;
occa::kernel gatherManyKernel_doubleAdd;
occa::kernel gatherManyKernel_doubleMul;
occa::kernel gatherManyKernel_doubleMin;
occa::kernel gatherManyKernel_doubleMax;
occa::kernel gatherManyKernel_intAdd;
occa::kernel gatherManyKernel_intMul;
occa::kernel gatherManyKernel_intMin;
occa::kernel gatherManyKernel_intMax;
occa::kernel gatherManyKernel_longAdd;
occa::kernel gatherManyKernel_longMul;
occa::kernel gatherManyKernel_longMin;
occa::kernel gatherManyKernel_longMax;

occa::kernel scatterKernel_float;
occa::kernel scatterKernel_double;
occa::kernel scatterKernel_int;
occa::kernel scatterKernel_long;
occa::kernel scatterVecKernel_float;
occa::kernel scatterVecKernel_double;
occa::kernel scatterVecKernel_int;
occa::kernel scatterVecKernel_long;
occa::kernel scatterManyKernel_float;
occa::kernel scatterManyKernel_double;
occa::kernel scatterManyKernel_int;
occa::kernel scatterManyKernel_long;
} // namespace ogs

void ogs::initKernels()
{

  auto buildKernel = [&](const std::string _fileName, const std::string& kernelName)
  {
    const auto& props = ogs::kernelInfo;

    const auto oklpath = std::string(getenv("NEKRS_KERNEL_DIR")) + "/core/ogs/";
    const auto fileName = oklpath + _fileName + ".okl";
    const auto reqName = fileName;

    if (platform->options.compareArgs("REGISTER ONLY", "TRUE")) {
      platform->kernelRequests.add(reqName, fileName, props); 
      return occa::kernel();
    } else {
      return platform->kernelRequests.load(reqName, kernelName);
    }
  };

  ogs::gatherScatterKernel_floatAdd =
      buildKernel("gatherScatter", "gatherScatter_float_add");
  ogs::gatherScatterNewKernel_floatAdd =
      buildKernel("gatherScatterNew", "gatherScatter_float_add");
  ogs::gatherScatterKernel_floatMul =
      buildKernel("gatherScatter", "gatherScatter_float_mul");
  ogs::gatherScatterKernel_floatMin =
      buildKernel("gatherScatter", "gatherScatter_float_min");
  ogs::gatherScatterKernel_floatMax =
      buildKernel("gatherScatter", "gatherScatter_float_max");

  ogs::gatherScatterKernel_doubleAdd =
      buildKernel("gatherScatter", "gatherScatter_double_add");
  ogs::gatherScatterNewKernel_doubleAdd =
      buildKernel("gatherScatterNew", "gatherScatter_double_add");
  ogs::gatherScatterKernel_doubleMul =
      buildKernel("gatherScatter", "gatherScatter_double_mul");
  ogs::gatherScatterKernel_doubleMin =
      buildKernel("gatherScatter", "gatherScatter_double_min");
  ogs::gatherScatterNewKernel_doubleMin =
      buildKernel("gatherScatterNew", "gatherScatter_double_min");
  ogs::gatherScatterKernel_doubleMax =
      buildKernel("gatherScatter", "gatherScatter_double_max");
  ogs::gatherScatterNewKernel_doubleMax =
      buildKernel("gatherScatterNew", "gatherScatter_double_max");

  ogs::gatherScatterKernel_intAdd = buildKernel("gatherScatter", "gatherScatter_int_add");
  ogs::gatherScatterKernel_intMul = buildKernel("gatherScatter", "gatherScatter_int_mul");
  ogs::gatherScatterKernel_intMin = buildKernel("gatherScatter", "gatherScatter_int_min");
  ogs::gatherScatterKernel_intMax = buildKernel("gatherScatter", "gatherScatter_int_max");

  ogs::gatherScatterKernel_longAdd =
      buildKernel("gatherScatter", "gatherScatter_long_long_add");
  ogs::gatherScatterKernel_longMul =
      buildKernel("gatherScatter", "gatherScatter_long_long_mul");
  ogs::gatherScatterKernel_longMin =
      buildKernel("gatherScatter", "gatherScatter_long_long_min");
  ogs::gatherScatterKernel_longMax =
      buildKernel("gatherScatter", "gatherScatter_long_long_max");

  // gatherScatterMany

  ogs::gatherScatterManyKernel_floatAdd =
      buildKernel("gatherScatterMany", "gatherScatterMany_float_add");
  ogs::gatherScatterManyKernel_doubleAdd =
      buildKernel("gatherScatterMany", "gatherScatterMany_double_add");
  ogs::gatherScatterManyKernel_intAdd =
      buildKernel("gatherScatterMany", "gatherScatterMany_int_add");
  ogs::gatherScatterManyKernel_longAdd =
      buildKernel("gatherScatterMany", "gatherScatterMany_long_long_add");


  ogs::gatherScatterManyKernel_floatMul =
      buildKernel("gatherScatterMany", "gatherScatterMany_float_mul");
  ogs::gatherScatterManyKernel_floatMin =
      buildKernel("gatherScatterMany", "gatherScatterMany_float_min");
  ogs::gatherScatterManyKernel_floatMax =
      buildKernel("gatherScatterMany", "gatherScatterMany_float_max");

  ogs::gatherScatterManyKernel_doubleMul =
      buildKernel("gatherScatterMany", "gatherScatterMany_double_mul");
  ogs::gatherScatterManyKernel_doubleMin =
      buildKernel("gatherScatterMany", "gatherScatterMany_double_min");
  ogs::gatherScatterManyKernel_doubleMax =
      buildKernel("gatherScatterMany", "gatherScatterMany_double_max");

  ogs::gatherScatterManyKernel_intMul =
      buildKernel("gatherScatterMany", "gatherScatterMany_int_mul");
  ogs::gatherScatterManyKernel_intMin =
      buildKernel("gatherScatterMany", "gatherScatterMany_int_min");
  ogs::gatherScatterManyKernel_intMax =
      buildKernel("gatherScatterMany", "gatherScatterMany_int_max");

  ogs::gatherScatterManyKernel_longMul =
      buildKernel("gatherScatterMany", "gatherScatterMany_long_long_mul");
  ogs::gatherScatterManyKernel_longMin =
      buildKernel("gatherScatterMany", "gatherScatterMany_long_long_min");
  ogs::gatherScatterManyKernel_longMax =
      buildKernel("gatherScatterMany", "gatherScatterMany_long_long_max");

  // gather

  ogs::gatherKernel_floatAdd = buildKernel("gather", "gather_float_add");
  ogs::gatherKernel_doubleAdd = buildKernel("gather", "gather_double_add");
  ogs::gatherKernel_intAdd = buildKernel("gather", "gather_int_add");
  ogs::gatherKernel_longAdd = buildKernel("gather", "gather_long_long_add");

  ogs::gatherKernel_floatMul = buildKernel("gather", "gather_float_mul");
  ogs::gatherKernel_floatMin = buildKernel("gather", "gather_float_min");
  ogs::gatherKernel_floatMax = buildKernel("gather", "gather_float_max");

  ogs::gatherKernel_doubleMul = buildKernel("gather", "gather_double_mul");
  ogs::gatherKernel_doubleMin = buildKernel("gather", "gather_double_min");
  ogs::gatherKernel_doubleMax = buildKernel("gather", "gather_double_max");

  ogs::gatherKernel_intMul = buildKernel("gather", "gather_int_mul");
  ogs::gatherKernel_intMin = buildKernel("gather", "gather_int_min");
  ogs::gatherKernel_intMax = buildKernel("gather", "gather_int_max");

  ogs::gatherKernel_longMul = buildKernel("gather", "gather_long_long_mul");
  ogs::gatherKernel_longMin = buildKernel("gather", "gather_long_long_min");
  ogs::gatherKernel_longMax = buildKernel("gather", "gather_long_long_max");

  // gatherMany

  ogs::gatherManyKernel_floatAdd = buildKernel("gatherMany", "gatherMany_float_add");
  ogs::gatherManyKernel_doubleAdd = buildKernel("gatherMany", "gatherMany_double_add");
  ogs::gatherManyKernel_intAdd = buildKernel("gatherMany", "gatherMany_int_add");
  ogs::gatherManyKernel_longAdd = buildKernel("gatherMany", "gatherMany_long_long_add");

  ogs::gatherManyKernel_floatMul = buildKernel("gatherMany", "gatherMany_float_mul");
  ogs::gatherManyKernel_floatMin = buildKernel("gatherMany", "gatherMany_float_min");
  ogs::gatherManyKernel_floatMax = buildKernel("gatherMany", "gatherMany_float_max");

  ogs::gatherManyKernel_doubleMul = buildKernel("gatherMany", "gatherMany_double_mul");
  ogs::gatherManyKernel_doubleMin = buildKernel("gatherMany", "gatherMany_double_min");
  ogs::gatherManyKernel_doubleMax = buildKernel("gatherMany", "gatherMany_double_max");

  ogs::gatherManyKernel_intMul = buildKernel("gatherMany", "gatherMany_int_mul");
  ogs::gatherManyKernel_intMin = buildKernel("gatherMany", "gatherMany_int_min");
  ogs::gatherManyKernel_intMax = buildKernel("gatherMany", "gatherMany_int_max");

  ogs::gatherManyKernel_longMul = buildKernel("gatherMany", "gatherMany_long_long_mul");
  ogs::gatherManyKernel_longMin = buildKernel("gatherMany", "gatherMany_long_long_min");
  ogs::gatherManyKernel_longMax = buildKernel("gatherMany", "gatherMany_long_long_max");

  // scatter

  ogs::scatterKernel_float = buildKernel("scatter", "scatter_float");
  ogs::scatterKernel_double = buildKernel("scatter", "scatter_double");
  ogs::scatterKernel_int = buildKernel("scatter", "scatter_int");
  ogs::scatterKernel_long = buildKernel("scatter", "scatter_long_long");

  // scatterMany

  ogs::scatterManyKernel_float = buildKernel("scatterMany", "scatterMany_float");
  ogs::scatterManyKernel_double = buildKernel("scatterMany", "scatterMany_double");
  ogs::scatterManyKernel_int = buildKernel("scatterMany", "scatterMany_int");
  ogs::scatterManyKernel_long = buildKernel("scatterMany", "scatterMany_long_long");
}

void ogs::freeKernels()
{

  ogs::gatherScatterKernel_floatAdd.free();
  ogs::gatherScatterKernel_floatMul.free();
  ogs::gatherScatterKernel_floatMin.free();
  ogs::gatherScatterKernel_floatMax.free();
  ogs::gatherScatterKernel_doubleAdd.free();
  ogs::gatherScatterKernel_doubleMul.free();
  ogs::gatherScatterKernel_doubleMin.free();
  ogs::gatherScatterKernel_doubleMax.free();
  ogs::gatherScatterKernel_intAdd.free();
  ogs::gatherScatterKernel_intMul.free();
  ogs::gatherScatterKernel_intMin.free();
  ogs::gatherScatterKernel_intMax.free();
  ogs::gatherScatterKernel_longAdd.free();
  ogs::gatherScatterKernel_longMul.free();
  ogs::gatherScatterKernel_longMin.free();
  ogs::gatherScatterKernel_longMax.free();
  ogs::gatherScatterVecKernel_floatAdd.free();
  ogs::gatherScatterVecKernel_floatMul.free();
  ogs::gatherScatterVecKernel_floatMin.free();
  ogs::gatherScatterVecKernel_floatMax.free();
  ogs::gatherScatterVecKernel_doubleAdd.free();
  ogs::gatherScatterVecKernel_doubleMul.free();
  ogs::gatherScatterVecKernel_doubleMin.free();
  ogs::gatherScatterVecKernel_doubleMax.free();
  ogs::gatherScatterVecKernel_intAdd.free();
  ogs::gatherScatterVecKernel_intMul.free();
  ogs::gatherScatterVecKernel_intMin.free();
  ogs::gatherScatterVecKernel_intMax.free();
  ogs::gatherScatterVecKernel_longAdd.free();
  ogs::gatherScatterVecKernel_longMul.free();
  ogs::gatherScatterVecKernel_longMin.free();
  ogs::gatherScatterVecKernel_longMax.free();
  ogs::gatherScatterManyKernel_floatAdd.free();
  ogs::gatherScatterManyKernel_floatMul.free();
  ogs::gatherScatterManyKernel_floatMin.free();
  ogs::gatherScatterManyKernel_floatMax.free();
  ogs::gatherScatterManyKernel_doubleAdd.free();
  ogs::gatherScatterManyKernel_doubleMul.free();
  ogs::gatherScatterManyKernel_doubleMin.free();
  ogs::gatherScatterManyKernel_doubleMax.free();
  ogs::gatherScatterManyKernel_intAdd.free();
  ogs::gatherScatterManyKernel_intMul.free();
  ogs::gatherScatterManyKernel_intMin.free();
  ogs::gatherScatterManyKernel_intMax.free();
  ogs::gatherScatterManyKernel_longAdd.free();
  ogs::gatherScatterManyKernel_longMul.free();
  ogs::gatherScatterManyKernel_longMin.free();
  ogs::gatherScatterManyKernel_longMax.free();

  ogs::gatherKernel_floatAdd.free();
  ogs::gatherKernel_floatMul.free();
  ogs::gatherKernel_floatMin.free();
  ogs::gatherKernel_floatMax.free();
  ogs::gatherKernel_doubleAdd.free();
  ogs::gatherKernel_doubleMul.free();
  ogs::gatherKernel_doubleMin.free();
  ogs::gatherKernel_doubleMax.free();
  ogs::gatherKernel_intAdd.free();
  ogs::gatherKernel_intMul.free();
  ogs::gatherKernel_intMin.free();
  ogs::gatherKernel_intMax.free();
  ogs::gatherKernel_longAdd.free();
  ogs::gatherKernel_longMul.free();
  ogs::gatherKernel_longMin.free();
  ogs::gatherKernel_longMax.free();
  ogs::gatherVecKernel_floatAdd.free();
  ogs::gatherVecKernel_floatMul.free();
  ogs::gatherVecKernel_floatMin.free();
  ogs::gatherVecKernel_floatMax.free();
  ogs::gatherVecKernel_doubleAdd.free();
  ogs::gatherVecKernel_doubleMul.free();
  ogs::gatherVecKernel_doubleMin.free();
  ogs::gatherVecKernel_doubleMax.free();
  ogs::gatherVecKernel_intAdd.free();
  ogs::gatherVecKernel_intMul.free();
  ogs::gatherVecKernel_intMin.free();
  ogs::gatherVecKernel_intMax.free();
  ogs::gatherVecKernel_longAdd.free();
  ogs::gatherVecKernel_longMul.free();
  ogs::gatherVecKernel_longMin.free();
  ogs::gatherVecKernel_longMax.free();
  ogs::gatherManyKernel_floatAdd.free();
  ogs::gatherManyKernel_floatMul.free();
  ogs::gatherManyKernel_floatMin.free();
  ogs::gatherManyKernel_floatMax.free();
  ogs::gatherManyKernel_doubleAdd.free();
  ogs::gatherManyKernel_doubleMul.free();
  ogs::gatherManyKernel_doubleMin.free();
  ogs::gatherManyKernel_doubleMax.free();
  ogs::gatherManyKernel_intAdd.free();
  ogs::gatherManyKernel_intMul.free();
  ogs::gatherManyKernel_intMin.free();
  ogs::gatherManyKernel_intMax.free();
  ogs::gatherManyKernel_longAdd.free();
  ogs::gatherManyKernel_longMul.free();
  ogs::gatherManyKernel_longMin.free();
  ogs::gatherManyKernel_longMax.free();

  ogs::scatterKernel_float.free();
  ogs::scatterKernel_double.free();
  ogs::scatterKernel_int.free();
  ogs::scatterKernel_long.free();
  ogs::scatterVecKernel_float.free();
  ogs::scatterVecKernel_double.free();
  ogs::scatterVecKernel_int.free();
  ogs::scatterVecKernel_long.free();
  ogs::scatterManyKernel_float.free();
  ogs::scatterManyKernel_double.free();
  ogs::scatterManyKernel_int.free();
  ogs::scatterManyKernel_long.free();

  ogs::o_haloBuf.free();
  ogs::haloBuf = NULL;
}
