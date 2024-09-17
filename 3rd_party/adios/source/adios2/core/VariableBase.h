/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * VariableBase.h Base class for Variable and VariableStruct types. Contains
 * common elements.
 *
 *  Created on: Feb 20, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_CORE_VARIABLEBASE_H_
#define ADIOS2_CORE_VARIABLEBASE_H_

/// \cond EXCLUDE_FROM_DOXYGEN
#include <float.h>
#include <iostream>
#include <limits.h>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>
/// \endcond

#include "adios2/common/ADIOSConfig.h"
#include "adios2/common/ADIOSTypes.h"
#include "adios2/core/Operator.h"

namespace adios2
{
namespace core
{

// forward declaration for reading streaming mode
class IO;
class Engine;

/** Base class for Variable<T> (primitives) and VariableStruct classes */
class VariableBase
{

public:
    /** unique identifier inside Method that creates a Variable */
    const std::string m_Name;

    /** primitive from <T> or compound from struct */
    const DataType m_Type;

    /** Variable -> sizeof(T),
     *  VariableStruct -> from constructor sizeof(struct) */
    const size_t m_ElementSize;

    /* User requested memory space */
    MemorySpace m_MemSpace = MemorySpace::Detect;
#if defined(ADIOS2_HAVE_KOKKOS) || defined(ADIOS2_HAVE_GPU_SUPPORT)
    ArrayOrdering m_BaseLayout;
    ArrayOrdering m_ArrayLayout = ArrayOrdering::Auto;
#endif

    ShapeID m_ShapeID = ShapeID::Unknown; ///< see shape types in ADIOSTypes.h
    size_t m_BlockID = 0;                 ///< current block ID for local variables, global = 0
    SelectionType m_SelectionType = SelectionType::BoundingBox;

    bool m_SingleValue = false; ///< true: single value, false: array
    Dims m_Shape;               ///< total dimensions across MPI
    Dims m_Start;               ///< starting point (offsets) in global shape
    Dims m_Count;               ///< dimensions from m_Start in global shape

    Dims m_MemoryStart; ///< start offset
    Dims m_MemoryCount; ///< local dimensions

    /** Global array was written as Joined array, so read accordingly */
    bool m_ReadAsJoined = false;

    /** Global array was written as Local value, so read accordingly */
    bool m_ReadAsLocalValue = false;

    /** For read mode, false: streaming */
    bool m_RandomAccess = true;

    /** used in streaming mode, true: first variable encounter, false: variable
     * already encountered in previous step */
    bool m_FirstStreamingStep = true;

    std::vector<std::shared_ptr<Operator>> m_Operations;

    size_t m_AvailableStepsStart = 0;
    size_t m_AvailableStepsCount = 0;

    size_t m_StepsStart = 0;
    size_t m_StepsCount = 1;

    size_t m_JoinedDimPos = 0; // the joined dimension in a  JoinedArray

    /** Index Metadata Position in a serial metadata buffer */
    size_t m_IndexStart = 0;

    Engine *m_Engine = nullptr;

    /** user requested accuracy */
    Accuracy m_AccuracyRequested = {0.0, 0.0, false};
    /** provided accuracy */
    Accuracy m_AccuracyProvided = {0.0, 0.0, false};

    /** Index to Step and blocks' (inside a step) characteristics position in a
     * serial metadata buffer
     * <pre>
     * key: step number (time_index in bp3 format)
     * value:  vector of block starts for that step
     * </pre>
     * */
    std::map<size_t, std::vector<size_t>> m_AvailableStepBlockIndexOffsets;

    std::map<size_t, Dims> m_AvailableShapes;

    std::set<std::string> m_PrefixedVariables;
    std::set<std::string> m_PrefixedAttributes;

    VariableBase(const std::string &name, const DataType type, const size_t elementSize,
                 const Dims &shape, const Dims &start, const Dims &count, const bool constantShape);

    virtual ~VariableBase() = default;

    /**
     * Returns the total number of elements
     * @return number of elements
     */
    size_t TotalSize() const noexcept;

#if defined(ADIOS2_HAVE_KOKKOS) || defined(ADIOS2_HAVE_GPU_SUPPORT)
    /**
     * Get the layout used by the user buffers
     * @return the layout used by the user buffers (RowMajor or ColumnMajor)
     */
    ArrayOrdering GetArrayLayout();

    /**
     * Set the layout used by the user buffers
     * @param the layout that will be used by future put/gets
     */
    void SetArrayLayout(const ArrayOrdering layout);
#endif

    /**
     * Get the memory space where a given buffers was allocated
     * @param pointer to the user data
     */
    MemorySpace GetMemorySpace(const void *ptr);

    /**
     * Set the memory space where user buffers will be allocated
     * @param the memory space where the expected buffers were allocated
     */
    void SetMemorySpace(const MemorySpace mem);

    /**
     * Set new shape
     * @param shape input shape to be applied to this variable
     */
    void SetShape(const adios2::Dims &shape);

    /**
     * Use at read only for local variables
     * @param blockID
     */
    void SetBlockSelection(const size_t blockID);

    /**
     * Set new start and count dimensions
     * @param boxDims = {start, count}
     */
    void SetSelection(const Box<Dims> &boxDims);

    /**
     * Set the steps for the variable. The pointer passed at
     * reading must be able to hold enough memory to store multiple steps in a
     * single read. For writing it changes the time step
     * @param boxSteps {startStep, countStep}
     */
    void SetStepSelection(const Box<size_t> &boxSteps);

    /**
     * Set local offset and dimensions to memory passed at Put
     */
    void SetMemorySelection(const Box<Dims> &memorySelection);

    /**
     * Sets the requested accuracy for the next read operation.
     * The actual accuracy after the read is provided in GetAccuracy()
     */
    void SetAccuracy(const adios2::Accuracy &a) noexcept;

    /**
     * Get the provided accuracy for the last read operation.
     * Most operations provide data as it was written, meaning that
     * error is reported as 0.0
     */
    adios2::Accuracy GetAccuracy() const noexcept;

    /** Return the requested accuracy set by user with SetAccuracy */
    adios2::Accuracy GetAccuracyRequested() const noexcept;

    size_t GetAvailableStepsStart() const;

    size_t GetAvailableStepsCount() const;

    /**
     * Adds an operation to this variable.
     * @param op reference to an Operator object
     * @param parameters operation specific parameters
     * @return operator handler
     */
    size_t AddOperation(std::shared_ptr<core::Operator> op) noexcept;

    size_t AddOperation(const std::string &op, const Params &parameters = Params()) noexcept;

    /**
     * Removes all current Operations associated with AddOperation.
     * Provides the posibility to apply or not operators on a step basis.
     */
    void RemoveOperations() noexcept;

    /**
     * Sets a parameter by key/value in an existing operation from AddOperation
     * @param operationID returned handler form AddOperation
     * @param key input parameter key
     * @param value input parameter value
     */
    void SetOperationParameter(const size_t operationID, const std::string key,
                               const std::string value);

    /** Self-check dims according to type, called from Engine before Write
     * @param hint extra debugging info for the exception */
    void CheckDimensions(const std::string hint) const;

    bool IsConstantDims() const noexcept;
    void SetConstantDims() noexcept;

    bool IsValidStep(const size_t step) const noexcept;

    /**
     * Resets m_StepsStart and m_StepsCount. Must be called in BeginStep
     */
    void ResetStepsSelection(const bool zeroStart) noexcept;

    /**
     * Checks if variable has a conflict to be accessed as a stream and
     * random-access (SetStepSelection has been called)
     * @param hint improve exception error message
     * @throws std::invalid_argument if random access and streaming are called
     */
    void CheckRandomAccessConflict(const std::string hint) const;

    Dims Shape(const size_t step, const MemorySpace memSpace,
               const ArrayOrdering layout = ArrayOrdering::Auto) const;
    Dims Shape(const size_t step = adios2::EngineCurrentStep) const;

    /**
     * Get info for attributes associated with this variable. Attribute name
     * must start with variable.m_Name + separator
     * @param io
     * @param separator
     * @return attributes info
     */
    std::map<std::string, Params> GetAttributesInfo(core::IO &io, const std::string separator,
                                                    const bool fullNameKeys) const noexcept;

protected:
    bool m_ConstantDims = false; ///< true: fix m_Shape, m_Start, m_Count

    unsigned int m_DeferredCounter = 0;

    void InitShapeType();

    /** Self-check dims according to type, called right after DefineVariable and
     *  SetSelection.
     * @param hint extra debugging info for the exception */
    void CheckDimensionsCommon(const std::string hint) const;

    void CheckRandomAccess(const size_t step, const std::string hint) const;

#if defined(ADIOS2_HAVE_KOKKOS) || defined(ADIOS2_HAVE_GPU_SUPPORT)
    inline void UpdateLayout(Dims &shape);
#endif
};

} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_CORE_VARIABLEBASE_H_ */
