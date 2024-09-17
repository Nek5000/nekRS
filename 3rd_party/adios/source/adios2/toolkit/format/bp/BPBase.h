/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BPBase.h
 *
 *  Created on: Sep 3, 2019
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_TOOLKIT_FORMAT_BP_BPBASE_H_
#define ADIOS2_TOOLKIT_FORMAT_BP_BPBASE_H_

#include <bitset>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include "adios2/common/ADIOSConfig.h"
#include "adios2/common/ADIOSMacros.h"
#include "adios2/common/ADIOSTypes.h"
#include "adios2/helper/adiosFunctions.h"
#include "adios2/toolkit/aggregator/mpi/MPIChain.h"
#include "adios2/toolkit/format/bp/bpBackCompatOperation/BPBackCompatOperation.h"
#include "adios2/toolkit/format/buffer/Buffer.h"
#include "adios2/toolkit/format/buffer/heap/BufferSTL.h"
#include "adios2/toolkit/profiling/iochrono/IOChrono.h"

namespace adios2
{
namespace format
{

/** Base class for BP3 and BP4 Serializers and Deserializers */
class BPBase
{
public:
    /**
     * Metadata index used for Variables and Attributes, needed in a
     * container for characteristic sets merge independently for each Variable
     * or Attribute
     */
    struct SerialElementIndex
    {
        /** buffer containing the metadata index, start with 500bytes */
        std::vector<char> Buffer;

        /** number of characteristics sets (time and spatial aggregation) */
        uint64_t Count = 0;

        /** unique ID assigned to each variable for counter */
        const uint32_t MemberID;

        /**
         * starting point for offsets characteristics update (used in
         * aggregation) */
        size_t LastUpdatedPosition = 0;

        /**
         * flag indicating whether the variable is valid or not.
         * if it's not valid (meaning the variable is not put at current step
         * even though it was put in previous steps),
         * the variable index will not be copied to the
         * metadata buffer when EndStep of CurrentStep is called
         */
        bool Valid = false;

        /**
         * Record the time step of previous block of each variable used to
         * decide whether we should create a new header of the variable index or
         * not. If it's a new step, we create a new header. Otherwise,
         * we don't */
        uint32_t CurrentStep = 0;

        size_t CurrentHeaderPosition = 0;

        /**
         * Default constructor holding 200 bytes of memory
         * @param memberID id input from bp format
         * @param bufferSize initial buffer size
         */
        SerialElementIndex(const uint32_t memberID, const size_t bufferSize = 200);
    };

    struct MetadataSet
    {
        /**
         * updated with EndStep, if append it will be updated to last,
         * starts with 1 in adios1 legacy
         */
        uint32_t TimeStep = 1;

        /** single buffer for PGIndex */
        SerialElementIndex PGIndex = SerialElementIndex(0);

        // no priority for now
        /** @brief key: variable name, value: bp metadata variable index */
        std::unordered_map<std::string, SerialElementIndex> VarsIndices;

        /** @brief key: attribute name, value: bp metadata attribute index */
        std::unordered_map<std::string, SerialElementIndex> AttributesIndices;

        /** Fixed size for mini footer, adding 28 bytes for ADIOS version */
        const unsigned int MiniFooterSize = 28 + 28;

        /** number of current PGs */
        uint64_t DataPGCount = 0;
        /** current PG initial ( relative ) position in data buffer */
        size_t DataPGLengthPosition = 0;
        /** number of variables in current PG */
        uint32_t DataPGVarsCount = 0;
        /** current PG variable count ( relative ) position */
        size_t DataPGVarsCountPosition = 0;
        /** true: currently writing to a pg, false: no current pg */
        bool DataPGIsOpen = false;

        /** Used at Read, steps start at zero */
        size_t StepsStart = 0;

        /** Used at Read, number of total steps */
        size_t StepsCount = 0;

        /** Similar to TimeStep, but uses uint64_t and start from zero. Used for
         * streaming a large number of steps */
        size_t CurrentStep = 0;

        /** position of pg index in current buffer */
        size_t PGIndexStart;

        /** position of var index in current buffer */
        size_t VariablesIndexStart;

        /** position of attr index in current buffer */
        size_t AttributesIndexStart;

        /** length of metadata file to which we append */
        size_t MetadataFileLength = 0;
    };

    struct Minifooter
    {
        std::string VersionTag;
        uint64_t PGIndexStart = 0;
        uint64_t VarsIndexStart = 0;
        uint64_t AttributesIndexStart = 0;
        int8_t Version = -1;
        uint8_t ADIOSVersionMajor = 0;
        uint8_t ADIOSVersionMinor = 0;
        uint8_t ADIOSVersionPatch = 0;
        uint32_t ADIOSVersion = 0; // major*1M + minor*1k + patch e.g. 2007001
        bool IsLittleEndian = true;
        bool HasSubFiles = false;

        Minifooter(const int8_t version);
    };

    /** groups all user-level parameters in a single struct */
    struct Parameters
    {
        /** Parameter to flush transports at every number of steps, to be used
         * at EndStep */
        size_t FlushStepsCount = 1;

        /** initial buffer size */
        size_t InitialBufferSize = DefaultInitialBufferSize;

        /** max buffer size */
        size_t MaxBufferSize = DefaultMaxBufferSize;

        /**
         * sub-block size for min/max calculation of large arrays in number of
         * elements (not bytes). The default big number per Put() default will
         * result in the original single min/max value-pair per block
         */
        size_t StatsBlockSize = 1125899906842624ULL;

        /** buffer memory growth factor */
        float GrowthFactor = DefaultBufferGrowthFactor;

        /** open timeout seconds BP4Only */
        float OpenTimeoutSecs = 0.f;

        /** Timeout for BeginStep in read. Set by the user as parameter
         * default value is calculated from the timeout value  */
        float BeginStepPollingFrequencySecs = 1.f;

        /** statistics verbosity, only 1 (BP4) is supported */
        unsigned int StatsLevel = 1;

        /** might be used in large payload copies to buffer */
        unsigned int Threads = 1;

        /** default time unit in m_Profiler */
        TimeUnit ProfileUnit = DefaultTimeUnitEnum;

        /** true: open files for write asynchronously,
         * false: all serial operations */
        bool AsyncOpen = true;

        /** true: write collective metadata, false: skip */
        bool CollectiveMetadata = true;

        /** true: NVMex each rank creates its own directory */
        bool NodeLocal = false;

        /** true: BeginStepPollingFrequency parameter is set */
        bool BeginStepPollingFrequencyIsSet = false;

        /** Burst buffer base path */
        std::string BurstBufferPath;

        /** Drain the file from Burst Buffer to the original path
         *  Relevant only if BurstBufferPath is set */
        bool BurstBufferDrain = true;
        /** Verbose level for burst buffer draining thread */
        int BurstBufferVerbose = 0;

        /** Stream reader flag: process metadata step-by-step
         * instead of parsing everything available
         */
        bool StreamReader = false;

        /** Number of aggregators.
         * Must be a value between 1 and number of MPI ranks
         * 0 as default means that the engine must define the number of
         * aggregators
         */
        unsigned int NumAggregators = 0;
    };

    /** Return type of the ResizeBuffer function. */
    enum class ResizeResult
    {
        Failure,   //!< FAILURE, caught a std::bad_alloc
        Unchanged, //!< UNCHANGED, no need to resize (sufficient capacity)
        Success,   //!< SUCCESS, resize was successful
        Flush      //!< FLUSH, need to flush to transports for current variable
    };

    helper::Comm const &m_Comm; ///< multi-process communicator from Engine
    int m_RankMPI = 0;          ///< current MPI rank process
    int m_SizeMPI = 1;          ///< current MPI processes size
    int m_Processes = 1;        ///< number of aggregated MPI processes

    /** used for appending to existing file */
    size_t m_PreMetadataFileLength = 0;
    /** used for appending to existing file */
    size_t m_PreDataFileLength = 0;

    /** contains data buffer for this rank */
    BufferSTL m_Data;

    /** contains collective metadata buffer, only used by rank 0 */
    BufferSTL m_Metadata;

    /** contains metadata indices */
    MetadataSet m_MetadataSet;

    /** contains user level parameters */
    Parameters m_Parameters;

    /** true: Close was called, Engine will call this many times for different
     * transports */
    bool m_IsClosed = false;

    /** buffering and MPI aggregation profiling info, set by user */
    profiling::IOChrono m_Profiler;

    /** from host language in data information at read */
    bool m_IsRowMajor = true;

    /** if reader and writer have different ordering (column vs row major) */
    bool m_ReverseDimensions = false;

    /** manages all communication tasks in aggregation */
    aggregator::MPIChain m_Aggregator;

    /** tracks Put and Get variables in deferred mode */
    std::set<std::string> m_DeferredVariables;

    /** tracks the overall size of deferred variables */
    size_t m_DeferredVariablesDataSize = 0;

    /**
     * attributes are serialized only once, this set contains the names of ones
     * already serialized.
     */
    std::unordered_set<std::string> m_SerializedAttributes;

    /**
     * scratch memory buffers used for management of temporary memory buffers
     * per thread.
     * This allows thread-safety mostly is deserialization.
     * Indices:
     * [threadID][bufferID]
     */
    std::map<size_t, std::map<size_t, std::vector<char>>> m_ThreadBuffers;

    /**
     * Default constructor
     * @param comm communicator from Engine
     */
    BPBase(helper::Comm const &comm);

    virtual ~BPBase() = default;

    /**
     * Init base don parameters passed from the user to IO
     * @param parameters input parameters
     */
    void Init(const Params &parameters, const std::string hint, const std::string engineType = "");

    /****************** NEED to check if some are virtual */

    /**
     * Resizes the data buffer to hold new dataIn size
     * @param dataIn input size for new data
     * @param hint extra messaging for exception handling
     * @return Failure, Unchanged, Success, Flush
     */
    ResizeResult ResizeBuffer(const size_t dataIn, const std::string hint);

    /**
     * Sets buffer's positions to zero and fill buffer with zero char
     * @param bufferSTL buffer to be reset
     * @param resetAbsolutePosition true: both bufferSTL.m_Position and
     * bufferSTL.m_AbsolutePosition set to 0,   false(default): only
     * bufferSTL.m_Position
     * is set to zero,
     */
    void ResetBuffer(Buffer &buffer, const bool resetAbsolutePosition = false,
                     const bool zeroInitialize = true);

    /** Delete buffer memory manually */
    void DeleteBuffers();

    size_t DebugGetDataBufferSize() const;

protected:
    /** file I/O method type, adios1 legacy, only POSIX and MPI_AGG are used */
    enum IO_METHOD
    {
        METHOD_UNKNOWN = -2,
        METHOD_NULL = -1,
        METHOD_MPI = 0,
        METHOD_DATATAP = 1, // OBSOLETE

        METHOD_POSIX = 2,
        METHOD_DATASPACES = 3,
        METHOD_VTK = 4, // non-existent

        METHOD_POSIX_ASCII = 5, // non-existent

        METHOD_MPI_CIO = 6, // OBSOLETE
        METHOD_PHDF5 = 7,
        METHOD_PROVENANCE = 8, // OBSOLETE
        METHOD_MPI_STRIPE = 9, // OBSOLETE
        METHOD_MPI_LUSTRE = 10,
        METHOD_MPI_STAGGER = 11, // OBSOLETE
        METHOD_MPI_AGG = 12,     // OBSOLETE
        METHOD_ADAPTIVE = 13,    // OBSOLETE
        METHOD_POSIX1 = 14,      // OBSOLETE
        METHOD_NC4 = 15,
        METHOD_MPI_AMR = 16,
        METHOD_MPI_AMR1 = 17, // OBSOLETE
        METHOD_FLEXPATH = 18,
        METHOD_NSSI_STAGING = 19,
        METHOD_NSSI_FILTER = 20,
        METHOD_DIMES = 21,
        METHOD_VAR_MERGE = 22,
        METHOD_MPI_BGQ = 23,
        METHOD_ICEE = 24,
        METHOD_COUNT = 25,
        METHOD_FSTREAM = 26,
        METHOD_FILE = 27,
        METHOD_ZMQ = 28,
        METHOD_MDTM = 29
    };

    /** DataTypes mapping in BP Format, adios1 legacy */
    enum DataTypes
    {
        type_unknown = -1, //!< type_unknown
        type_byte = 0,     //!< type_byte
        type_short = 1,    //!< type_short
        type_integer = 2,  //!< type_integer
        type_long = 4,     //!< type_long

        type_unsigned_byte = 50,    //!< type_unsigned_byte
        type_unsigned_short = 51,   //!< type_unsigned_short
        type_unsigned_integer = 52, //!< type_unsigned_integer
        type_unsigned_long = 54,    //!< type_unsigned_long

        type_real = 5,        //!< type_real or float
        type_double = 6,      //!< type_double
        type_long_double = 7, //!< type_long_double

        type_string = 9,          //!< type_string
        type_complex = 10,        //!< type_complex
        type_double_complex = 11, //!< type_double_complex
        type_string_array = 12,   //!< type_string_array

        type_char = 55 //!< type_char
    };

    /** Maps C++ type to DataTypes enum */
    template <typename T>
    struct TypeTraits;

    /** Characteristic ID in variable metadata, legacy adios1 */
    enum CharacteristicID
    {
        characteristic_value = 0,           //!< characteristic_value
        characteristic_min = 1,             //!< Used to read in older bp file format
        characteristic_max = 2,             //!< Used to read in older bp file format
        characteristic_offset = 3,          //!< characteristic_offset
        characteristic_dimensions = 4,      //!< characteristic_dimensions
        characteristic_var_id = 5,          //!< characteristic_var_id
        characteristic_payload_offset = 6,  //!< characteristic_payload_offset
        characteristic_file_index = 7,      //!< characteristic_file_index
        characteristic_time_index = 8,      //!< characteristic_time_index
        characteristic_bitmap = 9,          //!< characteristic_bitmap
        characteristic_stat = 10,           //!< characteristic_stat
        characteristic_transform_type = 11, //!< characteristic_transform_type
        characteristic_minmax = 12          //!< min-max array for subblocks
    };

    /** Define statistics type for characteristic ID = 10 */
    enum VariableStatistics
    {
        statistic_min = 0,
        statistic_max = 1,
        statistic_cnt = 2,
        statistic_sum = 3,
        statistic_sum_square = 4,
        statistic_hist = 5,
        statistic_finite = 6
    };

    /** Define transform types for characteristic ID = 11 */
    enum TransformTypes
    {
        transform_unknown = -1,
        transform_none = 0,
        transform_identity = 1,
        transform_zlib = 2,
        transform_bzip2 = 3,
        transform_szip = 4,
        transform_isobar = 5,
        transform_aplod = 6,
        transform_alacrity = 7,
        transform_zfp = 8,
        transform_sz = 9,
        transform_lz4 = 10,
        transform_blosc = 11,
        transform_mgard = 12,
        transform_png = 13,
        transform_sirius = 14,
        transform_mgardplus = 15,
        transform_plugin = 16,
    };

    /** Supported transform types */
    static const std::set<std::string> m_TransformTypes;

    /** Mapping of transform types enaum to string */
    static const std::map<int, std::string> m_TransformTypesToNames;

    // Transform related functions
    /**
     * Translates string to enum
     * @param transformType input
     * @return corresponding enum TransformTypes
     */
    TransformTypes TransformTypeEnum(const std::string transformType) const noexcept;

    /**
     * Returns the proper derived class for BPOperation based on type
     * @param type input, must be a supported type under BPOperation
     * @return derived class if supported, false pointer if type not supported
     */
    std::shared_ptr<BPBackCompatOperation>
    SetBPBackCompatOperation(const std::string type) const noexcept;

    struct ProcessGroupIndex
    {
        uint64_t Offset;
        uint32_t Step;
        int32_t ProcessID;
        uint16_t Length;
        std::string Name;
        std::string StepName;
        char IsColumnMajor;
    };

    /** holds extra metadata for operations in BP buffer */
    struct BPOpInfo
    {
        std::vector<char> Metadata;
        // pre-operator dimensions
        Dims PreShape;
        Dims PreCount;
        Dims PreStart;
        std::string Type; // Operator type, not data type
        uint8_t PreDataType;
        bool IsActive = false;
    };

    /** Fields for statistics characteristic ID = 10 */
    template <class T>
    struct Stats
    {
        std::vector<T> Values;
        std::vector<T> MinMaxs; // sub-block level min-max
        struct helper::BlockDivisionInfo SubBlockInfo;
        double BitSum = 0.;
        double BitSumSquare = 0.;
        uint64_t Offset = 0;
        uint64_t PayloadOffset = 0;
        T Min;
        T Max;
        T Value;
        uint32_t Step = 0;
        uint32_t FileIndex = 0;
        uint32_t MemberID = 0;
        uint32_t BitCount = 0;
        std::bitset<32> Bitmap;
        uint8_t BitFinite = 0;
        bool IsValue = false;
        BPOpInfo Op;

        Stats() : Min(), Max(), Value() {}
    };

    /**
     * Fields for all characteristics in BP buffer
     */
    template <class T>
    struct Characteristics
    {
        Stats<T> Statistics;
        Dims Shape;
        Dims Start;
        Dims Count;
        ShapeID EntryShapeID = ShapeID::Unknown;
        uint32_t EntryLength = 0;
        uint8_t EntryCount = 0;
    };

    struct ElementIndexHeader
    {
        uint64_t CharacteristicsSetsCount;
        uint32_t Length;
        uint32_t MemberID;
        std::string GroupName;
        std::string Name;
        std::string Path;
        uint8_t DataType = std::numeric_limits<uint8_t>::max() - 1;
        int8_t Order;
    };

    /**
     * Convert input transport types from user strings to enum IO_METHOD
     * @param transportsTypes input user transports
     * @return vector with enum IO_METHOD
     */
    std::vector<uint8_t>
    GetTransportIDs(const std::vector<std::string> &transportsTypes) const noexcept;

    /**
     * Calculates the Process Index size in bytes according to the BP
     * format,
     * including list of method with no parameters (for now)
     * @param name
     * @param timeStepName
     * @param transportsSize
     * @return size of pg index
     */
    size_t GetProcessGroupIndexSize(const std::string name, const std::string timeStepName,
                                    const size_t transportsSize) const noexcept;

    /**
     * Reads a PG index from a buffer position and advances the position until
     * done
     * @param buffer input buffer
     * @param position input start position, output as end PGIndex position
     * @param isLittleEndian true: buffer is little endian, false: big endian
     * @return populated PGIndex struct
     */
    ProcessGroupIndex ReadProcessGroupIndexHeader(const std::vector<char> &buffer, size_t &position,
                                                  const bool isLittleEndian = true) const noexcept;

    /**
     * Reads a PG index from a buffer position and advances the position until
     * done
     * @param buffer input buffer
     * @param position input start position, output as end PGIndex position
     * @param isLittleEndian true: buffer is little endian, false: big endian
     * @return populated PGIndex struct
     */
    virtual ElementIndexHeader
    ReadElementIndexHeader(const std::vector<char> &buffer, size_t &position,
                           const bool isLittleEndian = true) const noexcept = 0;

    /**
     * Read variable element (block) characteristics from a buffer position and
     * advancing the position until done
     * @param buffer input buffer
     * @param position input start position, output as end element
     * characteristic position
     * @param untilTimeStep, stop if time step characteristic is found
     * @return populated Characteristics<T> struct
     */
    template <class T>
    Characteristics<T> ReadElementIndexCharacteristics(const std::vector<char> &buffer,
                                                       size_t &position, const DataTypes dataType,
                                                       size_t &joinedArrayShapePos,
                                                       const bool untilTimeStep = false,
                                                       const bool isLittleEndian = true) const;

    /**
     * Common function to extract a BP standard string, 2 bytes for length +
     * contents from a buffer position advancing the position until done
     * @param buffer input buffer
     * @param position input start position, output as end element
     * @return populated string
     */
    std::string ReadBPString(const std::vector<char> &buffer, size_t &position,
                             const bool isLittleEndian = true) const noexcept;

private:
    /**
     * Specialized template that populates an existing characteristics
     * struct for a variable block
     */
    template <class T>
    void ParseCharacteristics(const std::vector<char> &buffer, size_t &position,
                              const DataTypes dataType, const bool untilTimeStep,
                              Characteristics<T> &characteristics, size_t &joinedArrayShapePos,
                              const bool isLittleEndian = true) const;
};

} // end namespace format
} // end namespace adios2

#include "BPBase.inl"

#endif /* ADIOS2_TOOLKIT_FORMAT_BP_BPBASE_H_ */
