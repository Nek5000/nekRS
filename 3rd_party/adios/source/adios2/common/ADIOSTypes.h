/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * ADIOSTypes.h : public header that contains "using/typedef" alias, defaults
 * and parameters options as enum classes
 *
 *  Created on: Mar 23, 2017
 *      Author: Chuck Atkins chuck.atkins@kitware.com
 *              Norbert Podhorszki pnorbert@ornl.gov
 *              William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_ADIOSTYPES_H_
#define ADIOS2_ADIOSTYPES_H_

/// \cond EXCLUDE_FROM_DOXYGEN
#include <cstddef>
#include <cstdint>

#include <complex>
#include <limits>
#include <map>
#include <string>
#include <type_traits>
#include <utility> //std::pair
#include <vector>
/// \endcond

#include "adios2/common/ADIOSConfig.h"

namespace adios2
{

#ifdef ADIOS2_HAVE_DERIVED_VARIABLE
/** Type of derived variables */
enum class DerivedVarType
{
    MetadataOnly,     ///< Store only the metadata (default)
    ExpressionString, ///< Store only the expression string
    StoreData         ///< Store data and metadata
};
#endif

/** Memory space for the user provided buffers */
enum class MemorySpace
{
    Detect, ///< Detect the memory space automatically
    Host,   ///< Host memory space
#ifdef ADIOS2_HAVE_GPU_SUPPORT
    GPU ///< GPU memory space
#endif
};

/** Variable shape type identifier, assigned automatically from the signature of
 *  DefineVariable */
enum class ShapeID
{
    Unknown,     ///< undefined shapeID
    GlobalValue, ///< single global value, common case
    GlobalArray, ///< global (across MPI_Comm) array, common case
    JoinedArray, ///< global array with a common (joinable) dimension
    LocalValue,  ///< special case, local independent value
    LocalArray,  ///< special case, local independent array
};

/** Used to set IO class */
enum class IOMode
{
    Independent, ///< all I/O operations are independent per rank
    Collective   ///< expect collective I/O operations
};

/** OpenMode in IO Open */
enum class Mode
{
    Undefined,
    // open modes
    Write,
    Read,
    Append,
    ReadRandomAccess, // reader random access mode
    // launch execution modes
    Sync,
    Deferred
};

enum class ReadMultiplexPattern
{
    GlobalReaders,
    RoundRobin,
    FirstInFirstOut,
    OpenAllSteps
};

enum class StreamOpenMode
{
    Wait,
    NoWait
};

enum class ReadMode
{
    NonBlocking,
    Blocking
};

enum class StepMode
{
    Append,
    Update, // writer advance mode
    Read,   // reader advance mode
};

enum class StepStatus
{
    OK,
    NotReady,
    EndOfStream,
    OtherError
};

enum class TimeUnit
{
    Microseconds,
    Milliseconds,
    Seconds,
    Minutes,
    Hours
};

/** Type of selection */
enum class SelectionType
{
    BoundingBox, ///< Contiguous block of data defined by offsets and counts
                 /// per dimension
    WriteBlock,  ///< Selection of an individual block written by a writer
                 /// process
};

// Data types.
enum class DataType
{
    None,
    Int8,
    Int16,
    Int32,
    Int64,
    UInt8,
    UInt16,
    UInt32,
    UInt64,
    Float,
    Double,
    LongDouble,
    FloatComplex,
    DoubleComplex,
    String,
    Char,
    Struct
};

/** Type of ArrayOrdering */
enum class ArrayOrdering
{
    RowMajor,    /// Contiguous elements of a row lie together in memory
    ColumnMajor, /// Contiguous elements of a column lie together in memory
    Auto         /// Default based on language type
};

// Types
using std::size_t;

using std::int8_t;
using std::int16_t;
using std::int32_t;
using std::int64_t;

using std::uint8_t;
using std::uint16_t;
using std::uint32_t;
using std::uint64_t;

// Complex
using cfloat = std::complex<float>;
using cdouble = std::complex<double>;

// Limit, using uint64_t to make it portable
constexpr uint64_t MaxU64 = std::numeric_limits<uint64_t>::max();
constexpr size_t MaxSizeT = std::numeric_limits<size_t>::max();
constexpr size_t DefaultSizeT = std::numeric_limits<size_t>::max();
constexpr size_t EngineCurrentStep = std::numeric_limits<size_t>::max();

union PrimitiveStdtypeUnion
{
    int8_t field_int8;
    uint8_t field_uint8;
    int16_t field_int16;
    uint16_t field_uint16;
    int32_t field_int32;
    uint32_t field_uint32;
    int64_t field_int64;
    uint64_t field_uint64;
    float field_float;
    double field_double;
    long double field_ldouble;
};

struct MinMaxStruct
{
    union PrimitiveStdtypeUnion MinUnion;
    union PrimitiveStdtypeUnion MaxUnion;
    void Init(DataType Type);
    void Dump(DataType Type);
};
struct MinBlockInfo
{
    int WriterID = 0;
    size_t BlockID = 0;
    size_t *Start;
    size_t *Count;
    MinMaxStruct MinMax;
    void *BufferP = NULL;
};

struct MinVarInfo
{
    size_t Step;
    bool WasLocalValue; // writer: localValue -> reader: 1D global array
    int Dims;
    const size_t *Shape;
    bool IsValue = false;
    bool IsReverseDims = false;
    std::vector<struct MinBlockInfo> BlocksInfo;
    MinVarInfo(const int D, const size_t *S)
    : Dims(D), Shape(S), IsValue(false), IsReverseDims(false), BlocksInfo({})
    {
    }
};

void PrintMVI(std::ostream &os, const MinVarInfo &mvi);

// adios defaults
const std::string DefaultFileLibrary("POSIX");
const std::string DefaultTimeUnit("Microseconds");
constexpr TimeUnit DefaultTimeUnitEnum(TimeUnit::Microseconds);

/** default initial bp buffer size, 16Kb, in bytes */
constexpr size_t DefaultInitialBufferSize = 16 * 1024;

/** default maximum bp buffer size, unlimited, in bytes.
 *  Needs to be studied for optimizing applications */
constexpr uint64_t DefaultMaxBufferSize = MaxSizeT - 1;

/** default buffer growth factor. Needs to be studied
 * for optimizing applications*/
constexpr float DefaultBufferGrowthFactor = 1.05f;

/** default Buffer Chunk Size
 *  128Mb */
constexpr uint64_t DefaultBufferChunkSize = 128 * 1024 * 1024;

/** default minimum size not copying deferred writes
 *  4Mb */
constexpr size_t DefaultMinDeferredSize = 4 * 1024 * 1024;

/** default size for writing/reading files using POSIX/fstream/stdio write
 *  2Gb - 100Kb (tolerance)*/
constexpr size_t DefaultMaxFileBatchSize = 2147381248;

/** default maximum shared memory segment size
 *  2 blocks of MaxFileBatchSize */
constexpr uint64_t DefaultMaxShmSize = 2 * DefaultMaxFileBatchSize;

constexpr char PathSeparator =
#ifdef _WIN32
    '\\';
#else
    '/';
#endif

// adios alias values and types
constexpr size_t UnknownDim = 0;
constexpr uint64_t JoinedDim = MaxU64 - 1;
constexpr uint64_t LocalValueDim = MaxU64 - 2;
constexpr uint64_t IrregularDim = MaxU64 - 3;
constexpr bool ConstantDims = true;
constexpr bool end_step = true;
constexpr bool LocalValue = true;
constexpr bool GlobalValue = false;

constexpr size_t UnknownStep = MaxU64;
constexpr double UnknownTime = std::numeric_limits<double>::infinity();

using Dims = std::vector<size_t>;
using Params = std::map<std::string, std::string>;
using vParams = std::vector<std::pair<std::string, Params>>;
using Steps = size_t;

template <class T>
using Box = std::pair<T, T>;

/**
 * TypeInfo
 * used to map from primitive types to stdint-based types
 */
template <typename T, typename Enable = void>
struct TypeInfo;

/** Data accuracy **/

/* Error. Accuracy can be requested for reading data.
   norm: 0.0 = L2, inf() = Linf
   relative: true = relative error, false = absolute error
 */
struct Accuracy
{
    double error;
    double norm;
    bool relative;
};

constexpr double L2_norm = 0.0;
constexpr double Linf_norm = std::numeric_limits<double>::infinity();

/**
 *  Return the actual size in bytes of elements of the given type.  Returns -1
 * for strings.
 */
int TypeElementSize(DataType adiosvartype);

/**
 *  Return true if a type has a reasonable Min/Max (is a simple integer or float type), false
 * otherwise
 */
bool TypeHasMinMax(DataType adiosvartype);

/**
 * ToString
 * makes a string from an enum class like ShapeID etc, for debugging etc
 * It is also overloaded elsewhere to allow for a readable representation of
 * Variable, Attribute, etc.
 */

std::string ToString(ShapeID value);
std::string ToString(IOMode value);
std::string ToString(Mode value);
std::string ToString(ReadMultiplexPattern value);
std::string ToString(StreamOpenMode value);
std::string ToString(ReadMode value);
std::string ToString(StepMode value);
std::string ToString(StepStatus value);
std::string ToString(TimeUnit value);
std::string ToString(SelectionType value);
std::string ToString(DataType type);
std::string ToString(const Dims &dims);
std::string ToString(const Box<Dims> &box);
std::string ToString(const MemorySpace value);

/** UserOptions holds all user options from ~/.config/adios2/adios2.yaml */
struct UserOptions
{
    struct General
    {
        int verbose;
    };

    struct Campaign
    {
        bool active;
        int verbose;
        std::string hostname;
        std::string campaignstorepath;
        std::string cachepath;
    };

    struct SST
    {
        int verbose;
    };

    General general;
    Campaign campaign;
    SST sst;
};

/**
 * os << [adios2_type] enables output of adios2 enums/classes directly
 * to output streams (e.g. std::cout), if ToString() can handle [adios2_type].
 */
template <typename T, typename Enable = decltype(ToString(std::declval<T>()))>
std::ostream &operator<<(std::ostream &os, const T &value);

namespace ops
{

// SZ PARAMETERS
#ifdef ADIOS2_HAVE_SZ

constexpr char LossySZ[] = "sz";

namespace sz
{
namespace key
{
constexpr char accuracy[] = "accuracy";
}
}

#endif

// ZFP PARAMETERS
#ifdef ADIOS2_HAVE_ZFP

constexpr char LossyZFP[] = "zfp";

namespace zfp
{
namespace key
{
constexpr char accuracy[] = "accuracy";
constexpr char backend[] = "backend";
constexpr char rate[] = "rate";
constexpr char precision[] = "precision";
}

namespace value
{
#ifdef ADIOS2_HAVE_ZFP_CUDA
constexpr char backend_cuda[] = "cuda";
#endif
constexpr char backend_omp[] = "omp";
constexpr char backend_serial[] = "serial";
}
}
#endif

// MGARD PARAMETERS
#ifdef ADIOS2_HAVE_MGARD

constexpr char LossyMGARD[] = "mgard";

namespace mgard
{
namespace key
{
constexpr char mode[] = "mode";
constexpr char tolerance[] = "tolerance";
constexpr char accuracy[] = "accuracy";
constexpr char s[] = "s";
}
}
#endif

#ifdef ADIOS2_HAVE_MGARD_MDR

constexpr char MDR[] = "mdr";

namespace mdr
{
constexpr double DOUBLE_ROUNDING_ERROR_LIMIT = 5.0e-16;
constexpr double FLOAT_ROUNDING_ERROR_LIMIT = 3.0e-7;
namespace key
{}
}
#endif

#ifdef ADIOS2_HAVE_LIBPRESSIO
constexpr char LossyLIBPRESSIO[] = "libpressio";
namespace libpressio
{
namespace key
{}
}
#endif

// PNG PARAMETERS
#ifdef ADIOS2_HAVE_PNG

constexpr char LosslessPNG[] = "png";

namespace png
{

namespace key
{
constexpr char color_type[] = "color_type";
constexpr char bit_depth[] = "bit_depth";
constexpr char compression_level[] = "compression_level";
}

namespace value
{
constexpr char color_type_GRAY[] = "PNG_COLOR_TYPE_GRAY";
constexpr char color_type_GRAY_ALPHA[] = "PNG_COLOR_TYPE_GRAY_ALPHA";
constexpr char color_type_PALETTE[] = "PNG_COLOR_TYPE_PALETTE";
constexpr char color_type_RGB[] = "PNG_COLOR_TYPE_RGB";
constexpr char color_type_RGB_ALPHA[] = "PNG_COLOR_TYPE_RGB_ALPHA";
constexpr char color_type_RGBA[] = "PNG_COLOR_TYPE_RGBA";
constexpr char color_type_GA[] = "PNG_COLOR_TYPE_GA";

constexpr char bit_depth_1[] = "1";
constexpr char bit_depth_2[] = "2";
constexpr char bit_depth_3[] = "3";
constexpr char bit_depth_4[] = "4";
constexpr char bit_depth_5[] = "5";
constexpr char bit_depth_6[] = "6";
constexpr char bit_depth_7[] = "7";
constexpr char bit_depth_8[] = "8";
constexpr char bit_depth_9[] = "9";

constexpr char compression_level_1[] = "1";
constexpr char compression_level_2[] = "2";
constexpr char compression_level_3[] = "3";
constexpr char compression_level_4[] = "4";
constexpr char compression_level_5[] = "5";
constexpr char compression_level_6[] = "6";
constexpr char compression_level_7[] = "7";
constexpr char compression_level_8[] = "8";
constexpr char compression_level_9[] = "9";

} // end namespace value

} // end namespace png
#endif

// BZIP2 PARAMETERS
#ifdef ADIOS2_HAVE_BZIP2

constexpr char LosslessBZIP2[] = "bzip2";
namespace bzip2
{

namespace key
{
constexpr char blockSize100k[] = "blockSize100k";
}

namespace value
{
constexpr char blockSize100k_1[] = "1";
constexpr char blockSize100k_2[] = "2";
constexpr char blockSize100k_3[] = "3";
constexpr char blockSize100k_4[] = "4";
constexpr char blockSize100k_5[] = "5";
constexpr char blockSize100k_6[] = "6";
constexpr char blockSize100k_7[] = "7";
constexpr char blockSize100k_8[] = "8";
constexpr char blockSize100k_9[] = "9";
} // end namespace value

} // end namespace bzip2
#endif

// Blosc PARAMETERS
#ifdef ADIOS2_HAVE_BLOSC2

constexpr char LosslessBlosc[] = "blosc";
namespace blosc
{

namespace key
{
constexpr char nthreads[] = "nthreads";
constexpr char compressor[] = "compressor";
constexpr char clevel[] = "clevel";
constexpr char doshuffle[] = "doshuffle";
constexpr char blocksize[] = "blocksize";
constexpr char threshold[] = "threshold";
}

namespace value
{

constexpr char compressor_blosclz[] = "blosclz";
constexpr char compressor_lz4[] = "lz4";
constexpr char compressor_lz4hc[] = "lz4hc";
constexpr char compressor_snappy[] = "snappy";
constexpr char compressor_zlib[] = "zlib";
constexpr char compressor_zstd[] = "zstd";

constexpr char clevel_0[] = "0";
constexpr char clevel_1[] = "1";
constexpr char clevel_2[] = "2";
constexpr char clevel_3[] = "3";
constexpr char clevel_4[] = "4";
constexpr char clevel_5[] = "5";
constexpr char clevel_6[] = "6";
constexpr char clevel_7[] = "7";
constexpr char clevel_8[] = "8";
constexpr char clevel_9[] = "9";

constexpr char doshuffle_shuffle[] = "BLOSC_SHUFFLE";
constexpr char doshuffle_noshuffle[] = "BLOSC_NOSHUFFLE";
constexpr char doshuffle_bitshuffle[] = "BLOSC_BITSHUFFLE";

} // end namespace value

} // end namespace blosc

#endif

} // end namespace ops

} // end namespace adios2

#include "ADIOSTypes.inl"

#endif /* ADIOS2_ADIOSTYPES_H_ */
