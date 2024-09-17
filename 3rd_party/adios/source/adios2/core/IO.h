/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IO.h factory class of Parameters, Variables, Transports to Engines
 *
 *  Created on: Dec 16, 2016
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_CORE_IO_H_
#define ADIOS2_CORE_IO_H_

/// \cond EXCLUDE_FROM_DOXYGEN
#include <map>
#include <memory> //std:shared_ptr
#include <string>
#include <unordered_map>
#include <utility> //std::pair
#include <vector>
/// \endcond

#include "adios2/common/ADIOSConfig.h"
#include "adios2/common/ADIOSMacros.h"
#include "adios2/common/ADIOSTypes.h"
#include "adios2/core/ADIOS.h"
#include "adios2/core/Attribute.h"
#include "adios2/core/CoreTypes.h"
#include "adios2/core/Group.h"
#include "adios2/core/Variable.h"
#ifdef ADIOS2_HAVE_DERIVED_VARIABLE
#include "adios2/core/VariableDerived.h"
#endif
#include "adios2/core/VariableStruct.h"

namespace adios2
{

namespace core
{

using VarMap = std::unordered_map<std::string, std::unique_ptr<VariableBase>>;
using AttrMap = std::unordered_map<std::string, std::unique_ptr<AttributeBase>>;

// forward declaration needed as IO is passed to Engine derived
// classes
class Engine;
class Group;

/** Factory class IO for settings, variables, and transports to an engine */
class IO
{

public:
    /** reference to object that created current IO */
    ADIOS &m_ADIOS;
    /** a pointer to a Group Object created from IO */
    std::shared_ptr<Group> m_Gr;

    /** unique identifier */
    const std::string m_Name;

    /** from ADIOS class passed to Engine created with Open */
    const std::string m_HostLanguage = "C++";

    ArrayOrdering m_ArrayOrder;

    /** From SetParameter, parameters for a particular engine from m_Type */
    Params m_Parameters;

    /** From AddTransport, parameters in map for each transport in vector */
    std::vector<Params> m_TransportsParameters;

    /** BP3 engine default if unknown */
    std::string m_EngineType = "File";

    /** at read for file engines: true: in streaming (step-by-step) mode, or
     * false: random-access mode (files) */
    bool m_ReadStreaming = false;

    /** used if m_Streaming is true by file reader engines */
    size_t m_EngineStep = 0;

    /** placeholder when reading XML file variable operations, executed until
     * DefineVariable in code */
    std::unordered_map<std::string, std::vector<std::pair<std::string, Params>>>
        m_VarOpsPlaceholder;

    /** true: prefix variables/attributes are cached per variable
     *   when function m_IsPrefixedNames called */
    bool m_IsPrefixedNames = false;

    /**
     * @brief Constructor called from ADIOS factory class DeclareIO function.
     * Not to be used direclty in applications.
     * @param adios reference to ADIOS object that owns current IO
     * @param name unique identifier for this IO object
     * @param inConfigFile IO defined in config file (XML)
     * @param hostLanguage current language using the adios2 library
     */
    IO(ADIOS &adios, const std::string name, const bool inConfigFile,
       const std::string hostLanguage);

    ~IO();

    IO(IO const &) = delete;
    IO &operator=(IO const &) = delete;

    /**
     * @brief Sets the engine type for this IO class object
     * @param engine predefined engine type, default is bpfile
     */
    void SetEngine(const std::string engine) noexcept;

    /**
     * @brief Set the IO mode (collective or independent), not yet implemented
     * @param IO mode */
    void SetIOMode(const IOMode mode);

    /**
     * @brief Version that passes a map to fill out parameters
     * initializer list = { "param1", "value1" },  {"param2", "value2"},
     * @param params adios::Params std::map<std::string, std::string>
     */
    void SetParameters(const Params &parameters = Params()) noexcept;

    /**
     * @brief Version that passes a single string to fill out many parameters.
     * initializer string = "param1=value1 , param2 = value2"
     * This function will throw std::invalid_argument for entries that
     * cannot be parsed into key=value pairs.
     */
    void SetParameters(const std::string &parameters);

    /**
     * @brief Sets a single parameter overwriting value if key exists;
     * @param key parameter key
     * @param value parameter value
     */
    void SetParameter(const std::string key, const std::string value) noexcept;

    /** @brief Retrieve current parameters map */
    Params &GetParameters() noexcept;

    /** @brief Delete all parameters */
    void ClearParameters() noexcept;

    /**
     * @brief Adds a transport and its parameters for the IO Engine
     * @param type must be a supported transport type
     * @param params acceptable parameters for a particular transport
     * @return transportIndex handler
     */
    size_t AddTransport(const std::string type, const Params &parameters = Params());

    /**
     * @brief Sets a single parameter to an existing transport identified with a
     * transportIndex handler from AddTransport.
     * This function overwrites existing parameter with the same key.
     * @param transportIndex index handler from AddTransport
     * @param key parameter key
     * @param value parameter value
     */
    void SetTransportParameter(const size_t transportIndex, const std::string key,
                               const std::string value);

    /**
     * @brief Define a Variable of primitive data type for current IO.
     * Default (name only) is a local single value,
     * in order to be compatible with ADIOS1.
     * @param name variable name, must be unique within Method
     * @param shape overall dimensions e.g. {Nx*size, Ny*size, Nz*size}
     * @param start point (offset) for MPI rank e.g. {Nx*rank, Ny*rank, Nz*rank}
     * @param count length for MPI rank e.g. {Nx, Ny, Nz}
     * @param constantShape true if dimensions, offsets and local sizes don't
     * change over time
     * @return reference to Variable object
     * @exception std::invalid_argument if Variable with unique name is already
     * defined
     */
    template <class T>
    Variable<T> &DefineVariable(const std::string &name, const Dims &shape = Dims(),
                                const Dims &start = Dims(), const Dims &count = Dims(),
                                const bool constantDims = false);
#ifdef ADIOS2_HAVE_DERIVED_VARIABLE
    VariableDerived &
    DefineDerivedVariable(const std::string &name, const std::string &expression,
                          const DerivedVarType varType = DerivedVarType::MetadataOnly);
#endif
    VariableStruct &DefineStructVariable(const std::string &name, StructDefinition &def,
                                         const Dims &shape = Dims(), const Dims &start = Dims(),
                                         const Dims &count = Dims(),
                                         const bool constantDims = false);

    /**
     * @brief Define array attribute
     * @param name must be unique for the IO object
     * @param array pointer to user data
     * @param elements number of data elements
     * @param variableName optionally associates the attribute to a Variable
     * @param allowModification true allows redefining/modifying existing
     * attribute
     * @return reference to internal Attribute
     * @exception std::invalid_argument if Attribute with unique name is already
     * defined
     */
    template <class T>
    Attribute<T> &DefineAttribute(const std::string &name, const T *array, const size_t elements,
                                  const std::string &variableName = "",
                                  const std::string separator = "/",
                                  const bool allowModification = false);

    /**
     * @brief Define single value attribute
     * @param name must be unique for the IO object
     * @param value single data value
     * @param allowModification true allows redefining/modifying existing
     * attribute
     * @return reference to internal Attribute
     * @exception std::invalid_argument if Attribute with unique name is already
     * defined
     */
    template <class T>
    Attribute<T> &
    DefineAttribute(const std::string &name, const T &value, const std::string &variableName = "",
                    const std::string separator = "/", const bool allowModification = false);

    /**
     * @brief Removes an existing Variable in current IO object.
     * Dangerous function since references and
     * pointers can be dangling after this call.
     * @param name unique identifier input
     * @return true: found and removed variable, false: not found, nothing to
     * remove
     */
    bool RemoveVariable(const std::string &name) noexcept;

    /**
     * @brief Removes all existing variables in current IO object.
     * Dangerous function since references and
     * pointers can be dangling after this call.
     */
    void RemoveAllVariables() noexcept;

    /**
     * @brief Removes an existing Attribute in current IO object.
     * Dangerous function since references and
     * pointers can be dangling after this call.
     * @param name unique identifier input
     * @return true: found and removed attribute, false: not found, nothing to
     * remove
     */
    bool RemoveAttribute(const std::string &name) noexcept;

    /**
     * @brief Removes all existing attributes in current IO object.
     * Dangerous function since references and
     * pointers can be dangling after this call.
     */
    void RemoveAllAttributes() noexcept;

    /**
     * @brief Retrieve map with variables info. Use when reading.
     * @param keys list of variable information keys to be extracted (case
     * insensitive). Leave empty to return all possible keys.
     * Possible values:
     * keys=['AvailableStepsCount','Type','Max','Min','SingleValue','Shape']
     * @return map with current variables and info
     * keys: Type, Min, Max, Value, AvailableStepsStart,
     * AvailableStepsCount, Shape, Start, Count, SingleValue
     */
    std::map<std::string, Params>
    GetAvailableVariables(const std::set<std::string> &keys = std::set<std::string>()) noexcept;

    /**
     * @brief Gets an existing variable of primitive type by name
     * @param name of variable to be retrieved
     * @return pointer to an existing variable in current IO, nullptr if not
     * found
     */
    template <class T>
    Variable<T> *InquireVariable(const std::string &name) noexcept;

    VariableStruct *InquireStructVariable(const std::string &name) noexcept;

    VariableStruct *InquireStructVariable(const std::string &name, const StructDefinition &def,
                                          const bool allowReorganize = false) noexcept;

    StructDefinition &DefineStruct(const std::string &name, const size_t size);

    /**
     * @brief Returns the type of an existing variable as an string
     * @param name input variable name
     * @return type primitive type
     */
    DataType InquireVariableType(const std::string &name) const noexcept;

    /**
     * Overload that accepts a const iterator into the m_Variables map if found
     * @param itVariable
     * @return type primitive type
     */
    DataType InquireVariableType(const VarMap::const_iterator itVariable) const noexcept;

    /**
     * Retrieves hash holding internal variable identifiers
     * @return
     * <pre>
     * key: unique variable name,
     * value: pointer to VariableBase
     * </pre>
     */
    const VarMap &GetVariables() const noexcept;
#ifdef ADIOS2_HAVE_DERIVED_VARIABLE
    const VarMap &GetDerivedVariables() const noexcept;
#endif

    /**
     * Retrieves hash holding internal Attributes identifiers
     * @return
     * <pre>
     * key: unique attribute name,
     * value: pointer to AttributeBase
     * </pre>
     */
    const AttrMap &GetAttributes() const noexcept;

    /**
     * Gets an existing attribute of primitive type by name
     * @param name of attribute to be retrieved
     * @return pointer to an existing attribute in current IO, nullptr if not
     * found
     */
    template <class T>
    Attribute<T> *InquireAttribute(const std::string &name, const std::string &variableName = "",
                                   const std::string separator = "/") noexcept;

    /**
     * @brief Returns the type of an existing attribute as an string
     * @param name input attribute name
     * @return type if found returns type as string, otherwise an empty string
     */
    DataType InquireAttributeType(const std::string &name, const std::string &variableName = "",
                                  const std::string separator = "/") const noexcept;

    /**
     * @brief Retrieve map with attributes info. Use when reading.
     * @return map with current attributes and info
     * keys: Type, Elements, Value
     */
    std::map<std::string, Params>
    GetAvailableAttributes(const std::string &variableName = std::string(),
                           const std::string separator = "/",
                           const bool fullNamesKeys = false) noexcept;

    /**
     * @brief Check existence in config file passed to ADIOS class constructor
     * @return true: defined in config file, false: not found in config file
     */
    bool InConfigFile() const noexcept;

    /**
     * Sets declared to true if IO exists in code created with ADIOS DeclareIO
     */
    void SetDeclared() noexcept;

    void SetArrayOrder(const ArrayOrdering ArrayOrder) noexcept;

    /**
     * Check if declared in code
     * @return true: created with ADIOS DeclareIO, false: dummy from config file
     */
    bool IsDeclared() const noexcept;

    /**
     * Adds an operator defined by the ADIOS class. Could be a variable set
     * transform, callback function, etc.
     * @param variable
     * @param operatorType
     * @param parameters
     */
    void AddOperation(const std::string &variable, const std::string &operatorType,
                      const Params &parameters = Params()) noexcept;

    /**
     * @brief Creates a polymorphic object that derives the Engine class,
     * based on the SetEngine function or config file input
     * @param name unique engine identifier within IO object
     * (e.g. file name in case of Files)
     * @param mode write, read, append from ADIOSTypes.h Mode
     * @param mpiComm assigns a new communicator to the Engine
     * @return a reference to a derived object of the Engine class
     * @exception std::invalid_argument if Engine with unique name is already
     * created with another Open
     */
    Engine &Open(const std::string &name, const Mode mode, helper::Comm comm);

    /**
     * Overloaded version that reuses the MPI_Comm object passed
     * from the ADIOS class to the IO class
     * @param name unique engine identifier within IO object
     * (file name in case of File transports)
     * @param mode write, read, append from ADIOSTypes.h OpenMode
     * @return a reference to a derived object of the Engine class
     * @exception std::invalid_argument if Engine with unique name is already
     * created with another Open
     */
    Engine &Open(const std::string &name, const Mode mode);

    /**
     * Retrieve an engine by name
     */
    Engine &GetEngine(const std::string &name);

    /**
     * Called from bindings Close function, not exposed in APIs.
     * Thin wrapper around map erase. Expected to call destructor on Engine in
     * map value.
     * @param name Engine name to be erased, must be created with Open
     */
    void RemoveEngine(const std::string &name);

    /**
     * Flushes all engines created with the current IO object using Open.
     * If no engine is created it does nothing.
     * @exception std::runtime_error if any engine Flush fails
     */
    void FlushAll();

    Group &CreateGroup(char delimiter);

    // READ FUNCTIONS, not yet implemented:
    /**
     * not yet implented
     * @param pattern
     */
    void SetReadMultiplexPattern(const ReadMultiplexPattern pattern);

    /**
     * not yet implemented
     * @param mode
     */
    void SetStreamOpenMode(const StreamOpenMode mode);

    /**
     * Resets all variables m_StepsStart and m_StepsCount
     * @param alwaysZero true: always m_StepsStart = 0, false: capture
     */
    void ResetVariablesStepSelection(const bool zeroStart = false, const std::string hint = "");

    void SetPrefixedNames(const bool isStep) noexcept;

    using MakeEngineFunc =
        std::function<std::shared_ptr<Engine>(IO &, const std::string &, const Mode, helper::Comm)>;
    struct EngineFactoryEntry
    {
        MakeEngineFunc MakeReader;
        MakeEngineFunc MakeWriter;
    };

    /**
     * Create a MakeEngineFunc that throws a std::invalid_argument
     * exception with the given error string.  This is useful when
     * an engine lacks support for either reading or writing.
     */
    static MakeEngineFunc NoEngine(std::string e);

    /**
     * Create an EngineFactoryEntry that throws a std::invalid_argument
     * exception with the given error string for both reader and writer.
     * This is useful when compiling without support for some engines.
     */
    static EngineFactoryEntry NoEngineEntry(std::string e);

    /**
     * Create an engine of type T.  This is intended to be used when
     * creating instances of EngineFactoryEntry for RegisterEngine.
     */
    template <typename T>
    static std::shared_ptr<Engine> MakeEngine(IO &io, const std::string &name, const Mode mode,
                                              helper::Comm comm)
    {
        return std::make_shared<T>(io, name, mode, std::move(comm));
    }

    /**
     * Register an engine factory entry to create a reader or writer
     * for an engine of the given engine type (named in lower case).
     */
    static void RegisterEngine(const std::string &engineType, EngineFactoryEntry entry);

    /*
     * Return list of all engines associated with this IO.
     */

    const std::map<std::string, std::shared_ptr<Engine>> &GetEngines() const { return m_Engines; }

    /** Inform about computation block through User->ADIOS */
    void EnterComputationBlock() noexcept;
    /** Inform about computation block through User->ADIOS */
    void ExitComputationBlock() noexcept;

private:
    /** true: exist in config file (XML) */
    const bool m_InConfigFile = false;

    bool m_IsDeclared = false;

    /** Independent (default) or Collective */
    adios2::IOMode m_IOMode = adios2::IOMode::Independent;

    VarMap m_Variables;
#ifdef ADIOS2_HAVE_DERIVED_VARIABLE
    VarMap m_VariablesDerived;
#endif

    AttrMap m_Attributes;

    std::map<std::string, std::shared_ptr<Engine>> m_Engines;

    /** Checks if attribute exists, called from DefineAttribute different
     *  signatures */
    void CheckAttributeCommon(const std::string &name) const;

    void CheckTransportType(const std::string type) const;

    template <class T>
    bool IsAvailableStep(const size_t step, const unsigned int variableIndex) noexcept;

    template <class T>
    Params GetVariableInfo(const std::string &variableName, const std::set<std::string> &keys);
};

} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_CORE_IO_H_ */
