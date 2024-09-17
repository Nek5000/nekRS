/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * VariableNT.h :
 *
 *  Created on: Apr 18, 2022
 *      Author: Jason Wang jason.ruonan.wang@gmail.com
 */

#ifndef ADIOS2_BINDINGS_CXX11_CXX11_VARIABLENT_H_
#define ADIOS2_BINDINGS_CXX11_CXX11_VARIABLENT_H_

#include "Operator.h"
#include "adios2/common/ADIOSTypes.h"
#include <string>

namespace adios2
{

namespace core
{
class VariableBase;
class StructDefinition;
}

class StructDefinition
{

public:
    void AddField(const std::string &name, const size_t offset, const DataType type,
                  const size_t size = 1);
    void Freeze() noexcept;
    size_t StructSize() const noexcept;
    std::string StructName() const noexcept;
    size_t Fields() const noexcept;
    std::string Name(const size_t index) const;
    size_t Offset(const size_t index) const;
    DataType Type(const size_t index) const;
    size_t ElementCount(const size_t index) const;
    explicit operator bool() const { return (m_StructDefinition != nullptr); }

private:
    friend class ADIOS;
    friend class IO;
    friend class VariableNT;
    StructDefinition(core::StructDefinition *ptr = nullptr);
    core::StructDefinition *m_StructDefinition = nullptr;
};

class VariableNT
{

public:
    VariableNT() = default;

    ~VariableNT() = default;

    /** Checks if object is valid, e.g. if( variable ) { //..valid } */
    explicit operator bool() const noexcept;

    /**
     * Sets the memory space for all following Puts
     * to either host (default) or device (currently only CUDA supported)
     * @param mem memory space where Put buffers are allocated
     */
    void SetMemorySpace(const MemorySpace mem);

    /**
     * Set new shape, care must be taken when reading back the variable for
     * different steps. Only applies to Global arrays.
     * @param shape new shape dimensions array
     */
    void SetShape(const adios2::Dims &shape);

    /**
     * Read mode only. Required for reading local variables, ShapeID() =
     * ShapeID::LocalArray or ShapeID::LocalValue. For Global Arrays it will Set
     * the appropriate Start and Count Selection for the global array
     * coordinates.
     * @param blockID: variable block index defined at write time. Blocks can be
     * inspected with bpls -D variableName
     */
    void SetBlockSelection(const size_t blockID);

    /**
     * Sets a variable selection modifying current {start, count}
     * Count is the dimension from Start point
     * @param selection input {start, count}
     */
    void SetSelection(const adios2::Box<adios2::Dims> &selection);

    /**
     * Set the local start (offset) point to the memory pointer passed at Put
     * and the memory local dimensions (count). Used for non-contiguous memory
     * writes and reads (e.g. multidimensional ghost-cells).
     * Currently Get only works for formats based on BP3.
     * @param memorySelection {memoryStart, memoryCount}
     * <pre>
     * 		memoryStart: relative local offset of variable.start to the
     * contiguous memory pointer passed at Put from which data starts. e.g. if
     * variable.Start() = {rank*Ny,0} and there is 1 ghost cell per dimension,
     * then memoryStart = {1,1}
     * 		memoryCount: local dimensions for the contiguous memory pointer
     * passed at Put, e.g. if there is 1 ghost cell per dimension and
     * variable.Count() = {Ny,Nx}, then memoryCount = {Ny+2,Nx+2}
     * </pre>
     */
    void SetMemorySelection(const adios2::Box<adios2::Dims> &memorySelection);

    /**
     * Sets a step selection modifying current startStep, countStep
     * countStep is the number of steps from startStep point
     * @param stepSelection input {startStep, countStep}
     */
    void SetStepSelection(const adios2::Box<size_t> &stepSelection);

    /**
     * Returns the number of elements required for pre-allocation based on
     * current count and stepsCount
     * @return elements of type T required for pre-allocation
     */
    size_t SelectionSize() const;

    /**
     * Inspects Variable name
     * @return name
     */
    std::string Name() const;

    /**
     * Inspects Variable type
     * @return type string literal containing the type: double, float, unsigned
     * int, etc.
     */
    std::string Type() const;

    /**
     * Inspects size of the current element type, sizeof(T)
     * @return sizeof(T) for current system
     */
    size_t Sizeof() const;

    /**
     * Inspects shape id for current variable
     * @return from enum adios2::ShapeID
     */
    adios2::ShapeID ShapeID() const;

    /**
     * Inspects shape in global variables
     * @param step input for a particular Shape if changing over time. If
     * default, either return absolute or in streaming mode it returns the shape
     * for the current engine step
     * @return shape vector
     */
    adios2::Dims Shape(const size_t step = adios2::EngineCurrentStep) const;

    /**
     * Inspects current start point
     * @return start vector
     */
    adios2::Dims Start() const;

    /**
     * Inspects current count from start
     * @return count vector
     */
    adios2::Dims Count() const;

    /**
     * For read mode, inspect the number of available steps
     * @return available steps
     */
    size_t Steps() const;

    /**
     * For read mode, inspect the start step for available steps
     * @return available start step
     */
    size_t StepsStart() const;

    /**
     * For read mode, retrieve current BlockID, default = 0 if not set with
     * SetBlockID
     * @return current block id
     */
    size_t BlockID() const;

    /**
     *Adds operation and parameters to current Variable object
     * @param op operator to be added
     * @param parameters key/value settings particular to the Variable, not to
     * be confused by op own parameters
     * @return operation index handler in Operations()
     */
    size_t AddOperation(const Operator op, const adios2::Params &parameters = adios2::Params());

    size_t AddOperation(const std::string &type,
                        const adios2::Params &parameters = adios2::Params());

    /**
     * Inspects current operators added with AddOperator
     * @return vector of Variable<T>::OperatorInfo
     */
    std::vector<Operator> Operations() const;

    /**
     * Removes all current Operations associated with AddOperation.
     * Provides the posibility to apply or not operators on a block basis.
     */
    void RemoveOperations();

    size_t StructFields() const;
    std::string StructFieldName(const size_t index) const;
    size_t StructFieldOffset(const size_t index) const;
    DataType StructFieldType(const size_t index) const;
    size_t StructFieldElementCount(const size_t index) const;

    union T
    {
        int8_t Int8;
        uint8_t UInt8;
        int16_t Int16;
        uint16_t UInt16;
        int32_t Int32;
        uint32_t UInt32;
        int64_t Int64;
        uint64_t UInt64;
        float Float;
        double Double;
        std::complex<float> Complex;
        std::complex<double> DComplex;
    };

    struct Info
    {
        adios2::Dims Start;
        adios2::Dims Count;
        int WriterID = 0;
        size_t BlockID = 0;
        size_t Step = 0;
        bool IsReverseDims = false;
        bool IsValue = false;
    };

    T Min(const size_t step = adios2::DefaultSizeT) const;
    T Max(const size_t step = adios2::DefaultSizeT) const;
    std::pair<T, T> MinMax(const size_t step = adios2::DefaultSizeT) const;

    double MinDouble(const size_t step = adios2::DefaultSizeT) const;
    double MaxDouble(const size_t step = adios2::DefaultSizeT) const;
    std::pair<double, double> MinMaxDouble(const size_t step = adios2::DefaultSizeT) const;

    StructDefinition GetWriteStructDef() noexcept;
    StructDefinition GetReadStructDef() noexcept;
    void SetReadStructDef(const StructDefinition &def);

private:
    friend class IO;
    friend class Engine;
    VariableNT(core::VariableBase *variable);
    core::VariableBase *m_Variable = nullptr;
};

std::string ToString(const VariableNT &variable);

} // end namespace adios2

#endif // ADIOS2_BINDINGS_CXX11_CXX11_VARIABLE_H_
