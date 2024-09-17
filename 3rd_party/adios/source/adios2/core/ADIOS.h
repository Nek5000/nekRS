/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * ADIOS.h : ADIOS library starting point, factory class for IO objects
 *  Created on: Oct 3, 2016
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_CORE_ADIOS_H_
#define ADIOS2_CORE_ADIOS_H_

/// \cond EXCLUDE_FROM_DOXYGEN
#include <functional> //std::function
#include <memory>     //std::shared_ptr
#include <string>
#include <unordered_map>
#include <vector>
/// \endcond

#include "adios2/common/ADIOSConfig.h"
#include "adios2/common/ADIOSTypes.h"
#include "adios2/core/CoreTypes.h"
#include "adios2/core/Operator.h"
#include "adios2/core/VariableStruct.h"
#include "adios2/helper/adiosComm.h"

// Campaign Manager as Global Service
#include <adios2/engine/campaign/CampaignManager.h>

namespace adios2
{
namespace core
{
class IO;

/** @brief Point of entry class for an application.
 *         Serves as factory of IO class objects and Transforms */
class ADIOS
{
public:
    /** Get the communicator passed to constructor for parallel case.  */
    helper::Comm const &GetComm() const { return m_Comm; }

    /** Changed by language bindings in constructor */
    const std::string m_HostLanguage = "C++";

    /**
     * @brief Constructor for MPI applications WITH a XML config file
     * @param configFile XML format (maybe support different formats in the
     * future (json)?)
     * @param mpiComm MPI communicator from application, make sure is valid
     * through the scope of adios2 calls
     */
    ADIOS(const std::string configFile, helper::Comm comm, const std::string hostLanguage);

    /**
     * @brief Constructor for non-MPI applications WITH a XML config file (it
     * must end with extension .xml)
     * @param configFile XML format (maybe support different formats in the
     * future (json)?)
     */
    ADIOS(const std::string configFile, const std::string hostLanguage);

    /**
     * @brief Constructor for MPI apps WITHOUT a XML config file
     * @param mpiComm MPI communicator from application
     */
    ADIOS(helper::Comm comm, const std::string hostLanguage);

    /**
     *  @brief ADIOS no-MPI default empty constructor
     */
    ADIOS(const std::string hostLanguage);

    /**
     * Delete copy constructor explicitly. Objects shouldn't be allowed to be
     * redefined. Use smart pointers if this is absolutely necessary.
     * @param adios reference to another adios object
     */
    ADIOS(const ADIOS &adios) = delete;

    ~ADIOS();

    /**
     * Declares a new IO class object and returns a reference to that object.
     * @param ioName must be unique
     * @return reference to newly created IO object inside current ADIOS object
     * @exception std::invalid_argument if IO with unique name is already
     * declared
     */
    IO &DeclareIO(const std::string name, const ArrayOrdering ArrayOrder = ArrayOrdering::Auto);

    /**
     * Retrieve a reference to an existing IO object created with DeclareIO.
     * Follow the C++11 STL containers at function.
     * @param name of IO to look for
     * @return if IO exists returns a reference to existing IO object inside
     * ADIOS
     * @exception std::invalid_argument if IO was not created with DeclareIO
     */
    IO &AtIO(const std::string name);

    /**
     * Flushes all engines in all IOs created with the current ADIOS object
     * using DeclareIO and IO.Open.
     * If no IO or engine is created it does nothing.
     * @exception std::runtime_error if any engine Flush fails
     */
    void FlushAll();

    /**
     * Declares a derived class of the Operator abstract class. If object is
     * defined in the user
     * config file, by name, it will be already created during the processing of
     * the config file. So this function returns a reference to that object.
     * @param name must be unique for each operator created with DefineOperator
     * @param type from derived class
     * @param parameters optional parameters
     * @return reference to Operator object
     * @exception std::invalid_argument if Operator with unique name is already
     * defined
     */
    std::pair<std::string, Params> &DefineOperator(const std::string &name, const std::string type,
                                                   const Params &parameters = Params());
    /**
     * Retrieve a reference pointer to an existing Operator object
     * created with DefineOperator.
     * @return if IO exists returns a reference to existing IO object inside
     * ADIOS, otherwise a nullptr
     */
    std::pair<std::string, Params> *InquireOperator(const std::string &name) noexcept;

    /*
     * StructDefinitions are defined using the operators in the IO,
     * but they are emplaced into a set in the ADIOS to give them
     * global scope
     */
    std::unordered_multimap<std::string, StructDefinition> m_StructDefinitions;

    /**
     * DANGER ZONE: removes a particular IO. This will effectively eliminate any
     * parameter from the config.xml file
     * @param name io input name
     * @return true: IO was found and removed, false: IO not found and not
     * removed
     */
    bool RemoveIO(const std::string name);

    /**
     * DANGER ZONE: removes all IOs created with DeclareIO. This will
     * effectively eliminate any parameter from the config.xml file
     */
    void RemoveAllIOs() noexcept;

    /** Inform ADIOS about entering communication-free computation block
     * in main thread. Useful when using Async IO */
    void EnterComputationBlock() noexcept;

    /** Inform ADIOS about exiting communication-free computation block
     * in main thread. Useful when using Async IO */
    void ExitComputationBlock() noexcept;

    void RecordOutputStep(const std::string &name, const size_t step = UnknownStep,
                          const double time = UnknownTime);

    /** A constant reference to the user options from ~/.config/adios2/adios2.yaml */
    const adios2::UserOptions &GetUserOptions();

private:
    /** Communicator given to parallel constructor. */
    helper::Comm m_Comm;

    /** XML File to be read containing configuration information */
    const std::string m_ConfigFile;
    std::string m_ConfigFileContents;

    /**
     * @brief List of IO class objects defined from either ADIOS
     * configuration file (XML) or the DeclareIO function explicitly.
     * Using map (binary tree) to preserve references returned by DeclareIO.
     * <pre>
     *     Key: unique method name
     *     Value: IO class object
     * </pre>
     */
    std::map<std::string, IO> m_IOs;

    /** operators created with DefineOperator */
    std::unordered_map<std::string, std::pair<std::string, Params>> m_Operators;

    /** campaign manager */
    engine::CampaignManager m_CampaignManager;

    /** Flag for Enter/ExitComputationBlock */
    bool enteredComputationBlock = false;

    void CheckOperator(const std::string name) const;

    std::string XMLInit(const std::string &configFileXML);

    void XMLIOInit(const std::string &configFileXML, const std::string &configFileContents,
                   core::IO &io);

    std::string YAMLInit(const std::string &configFileYAML);

    void YAMLInitIO(const std::string &configFileYAML, const std::string &configFileContents,
                    core::IO &io);

    adios2::UserOptions m_UserOptions;
    void SetUserOptionDefaults();
    void ProcessUserConfig();

private:
    /* Global services that we want to initialize at most once and shutdown
       automatically when the ADIOS object is destructed. This only works
       properly if the app creates an ADIOS object that is created before all
       other ADIOS objects and is destructed after all other ADIOS objects are
       destructed*/
    class GlobalServices;
    static class GlobalServices m_GlobalServices;

public:
    /** Global service AWS SDK initialization */
    static void Global_init_AWS_API();
};

} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_ADIOS_H_ */
