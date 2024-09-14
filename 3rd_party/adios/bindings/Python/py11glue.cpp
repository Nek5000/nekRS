/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * py11glue.cpp
 *
 *  Created on: Mar 16, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <complex>
#include <sstream>
#include <stdexcept>

#include <adios2.h>

#if ADIOS2_USE_MPI
#include <mpi4py/mpi4py.h>
#endif

#include "py11ADIOS.h"
#include "py11Attribute.h"
#include "py11Engine.h"
#include "py11IO.h"
#include "py11Operator.h"
#include "py11Query.h"
#include "py11Variable.h"

#if ADIOS2_USE_MPI

namespace pybind11
{
namespace detail
{
template <>
struct type_caster<adios2::py11::MPI4PY_Comm>
{
public:
    /**
     * This macro establishes the name 'MPI4PY_Comm' in
     * function signatures and declares a local variable
     * 'value' of type MPI4PY_Comm
     */
    PYBIND11_TYPE_CASTER(adios2::py11::MPI4PY_Comm, _("MPI4PY_Comm"));

    /**
     * Conversion part 1 (Python->C++): convert a PyObject into a MPI4PY_Comm
     * instance or return false upon failure. The second argument
     * indicates whether implicit conversions should be applied.
     */
    bool load(handle src, bool)
    {
        // Import mpi4py if it does not exist.
        if (!PyMPIComm_Get)
        {
            if (import_mpi4py() < 0)
            {
                throw std::runtime_error("ERROR: mpi4py not loaded correctly\n"); /* Python 2.X */
            }
        }
        // If src is not actually a MPI4PY communicator, the next
        // call returns nullptr, and we return false to indicate the conversion
        // failed.
        MPI_Comm *mpiCommPtr = PyMPIComm_Get(src.ptr());
        if (mpiCommPtr == nullptr)
        {
            return false;
        }
        value.comm = *mpiCommPtr;
        return true;
    }
};
} // namespace detail
} // namespace pybind11

#endif

PYBIND11_MODULE(ADIOS2_PYTHON_MODULE_NAME, m)
{
    m.attr("ConstantDims") = true;
    m.attr("VariableDims") = false;
    m.attr("LocalValueDim") = adios2::LocalValueDim;
    m.attr("GlobalValue") = false;
    m.attr("LocalValue") = true;

    m.attr("__version__") = ADIOS2_VERSION_STR;
#if defined(ADIOS2_HAVE_MPI)
    m.attr("is_built_with_mpi") = true;
#else
    m.attr("is_built_with_mpi") = false;
#endif
    m.attr("is_char_signed") = (char)-1 < 0;

    // enum classes
    pybind11::enum_<adios2::Mode>(m, "Mode")
        .value("Write", adios2::Mode::Write)
        .value("Read", adios2::Mode::Read)
        .value("ReadRandomAccess", adios2::Mode::ReadRandomAccess)
        .value("Append", adios2::Mode::Append)
        .value("Deferred", adios2::Mode::Deferred)
        .value("Sync", adios2::Mode::Sync)
        .export_values();

    pybind11::enum_<adios2::ShapeID>(m, "ShapeID")
        .value("Unknown", adios2::ShapeID::Unknown)
        .value("GlobalValue", adios2::ShapeID::GlobalValue)
        .value("GlobalArray", adios2::ShapeID::GlobalArray)
        .value("LocalValue", adios2::ShapeID::LocalValue)
        .value("LocalArray", adios2::ShapeID::LocalArray)
        .export_values();

    pybind11::enum_<adios2::StepMode>(m, "StepMode")
        .value("Append", adios2::StepMode::Append)
        .value("Update", adios2::StepMode::Update)
        .value("Read", adios2::StepMode::Read)
        .export_values();

    pybind11::enum_<adios2::StepStatus>(m, "StepStatus")
        .value("OK", adios2::StepStatus::OK)
        .value("NotReady", adios2::StepStatus::NotReady)
        .value("EndOfStream", adios2::StepStatus::EndOfStream)
        .value("OtherError", adios2::StepStatus::OtherError)
        .export_values();

    pybind11::class_<adios2::py11::ADIOS>(m, "ADIOS")
        // Python 2
        .def("__nonzero__",
             [](const adios2::py11::ADIOS &adios) {
                 const bool opBool = adios ? true : false;
                 return opBool;
             })
        // Python 3
        .def("__bool__",
             [](const adios2::py11::ADIOS &adios) {
                 const bool opBool = adios ? true : false;
                 return opBool;
             })
        .def(pybind11::init(), "adios2 module starting point "
                               "non-MPI, constructs an ADIOS class "
                               "object")
        .def(pybind11::init<const std::string &>(),
             "adios2 module starting point non-MPI, constructs an ADIOS class "
             "object",
             pybind11::arg("configFile"))
#if ADIOS2_USE_MPI
        .def(pybind11::init<const adios2::py11::MPI4PY_Comm>(),
             "adios2 module starting point, constructs an ADIOS class object",
             pybind11::arg("comm"))
        .def(pybind11::init<const std::string &, const adios2::py11::MPI4PY_Comm>(),
             "adios2 module starting point, constructs an ADIOS class object",
             pybind11::arg("configFile"), pybind11::arg("comm"))
#endif
        .def("DeclareIO", &adios2::py11::ADIOS::DeclareIO,
             "spawn IO object component returning a IO object with a unique "
             "name, throws an exception if IO with the same name is declared "
             "twice")
        .def("AtIO", &adios2::py11::ADIOS::AtIO,
             "returns an IO object "
             "previously defined IO object "
             "with DeclareIO, throws "
             "an exception if not found")
        .def("DefineOperator", &adios2::py11::ADIOS::DefineOperator)
        .def("InquireOperator", &adios2::py11::ADIOS::InquireOperator)
        .def("FlushAll", &adios2::py11::ADIOS::FlushAll,
             "flushes all engines in all spawned IO objects")
        .def("RemoveIO", &adios2::py11::ADIOS::RemoveIO,
             "DANGER ZONE: remove a particular IO by name, creates dangling "
             "objects to parameters, variable, attributes, engines created "
             "with removed IO")
        .def("RemoveAllIOs", &adios2::py11::ADIOS::RemoveAllIOs,
             "DANGER ZONE: remove all IOs in current ADIOS object, creates "
             "dangling objects to parameters, variable, attributes, engines "
             "created with removed IO");

    pybind11::class_<adios2::py11::IO>(m, "IO")
        // Python 2
        .def("__nonzero__",
             [](const adios2::py11::IO &io) {
                 const bool opBool = io ? true : false;
                 return opBool;
             })
        // Python 3
        .def("__bool__",
             [](const adios2::py11::IO &io) {
                 const bool opBool = io ? true : false;
                 return opBool;
             })
        .def("SetEngine", &adios2::py11::IO::SetEngine)
        .def("SetParameters", &adios2::py11::IO::SetParameters,
             pybind11::arg("parameters") = adios2::Params())
        .def("SetParameter", &adios2::py11::IO::SetParameter)
        .def("Parameters", &adios2::py11::IO::Parameters)
        .def("AddTransport", &adios2::py11::IO::AddTransport, pybind11::arg("type"),
             pybind11::arg("parameters") = adios2::Params())

        .def("DefineVariable",
             (adios2::py11::Variable(adios2::py11::IO::*)(
                 const std::string &, const pybind11::array &, const adios2::Dims &,
                 const adios2::Dims &, const adios2::Dims &, const bool)) &
                 adios2::py11::IO::DefineVariable,
             pybind11::return_value_policy::move, pybind11::arg("name"), pybind11::arg("array"),
             pybind11::arg("shape") = adios2::Dims(), pybind11::arg("start") = adios2::Dims(),
             pybind11::arg("count") = adios2::Dims(), pybind11::arg("isConstantDims") = false)

        .def("DefineVariable",
             (adios2::py11::Variable(adios2::py11::IO::*)(
                 const std::string &, const pybind11::object &, const adios2::Dims &,
                 const adios2::Dims &, const adios2::Dims &, const bool)) &
                 adios2::py11::IO::DefineVariable,
             pybind11::return_value_policy::move, pybind11::arg("name"), pybind11::arg("value"),
             pybind11::arg("shape") = adios2::Dims(), pybind11::arg("start") = adios2::Dims(),
             pybind11::arg("count") = adios2::Dims(), pybind11::arg("isConstantDims") = false)

        .def("DefineVariable",
             (adios2::py11::Variable(adios2::py11::IO::*)(const std::string &)) &
                 adios2::py11::IO::DefineVariable,
             pybind11::return_value_policy::move, pybind11::arg("name"))

        .def("InquireVariable", &adios2::py11::IO::InquireVariable,
             pybind11::return_value_policy::move)

        .def("InquireAttribute",
             (adios2::py11::Attribute(adios2::py11::IO::*)(const std::string &, const std::string &,
                                                           const std::string)) &
                 adios2::py11::IO::InquireAttribute,
             pybind11::arg("name"), pybind11::arg("variable_name") = "",
             pybind11::arg("separator") = "/", pybind11::return_value_policy::move)

        .def("DefineAttribute",
             (adios2::py11::Attribute(adios2::py11::IO::*)(
                 const std::string &, const pybind11::array &, const std::string &,
                 const std::string)) &
                 adios2::py11::IO::DefineAttribute,
             pybind11::arg("name"), pybind11::arg("array"), pybind11::arg("variable_name") = "",
             pybind11::arg("separator") = "/", pybind11::return_value_policy::move)

        .def(
            "DefineAttribute",
            (adios2::py11::Attribute(adios2::py11::IO::*)(const std::string &, const std::string &,
                                                          const std::string &, const std::string)) &
                adios2::py11::IO::DefineAttribute,
            pybind11::arg("name"), pybind11::arg("stringValue"),
            pybind11::arg("variable_name") = "", pybind11::arg("separator") = "/",
            pybind11::return_value_policy::move)

        .def("DefineAttribute",
             (adios2::py11::Attribute(adios2::py11::IO::*)(
                 const std::string &, const std::vector<std::string> &, const std::string &,
                 const std::string)) &
                 adios2::py11::IO::DefineAttribute,
             pybind11::arg("name"), pybind11::arg("strings"), pybind11::arg("variable_name") = "",
             pybind11::arg("separator") = "/", pybind11::return_value_policy::move)

        .def("DefineAttribute",
             (adios2::py11::Attribute(adios2::py11::IO::*)(
                 const std::string &, const std::vector<int> &, const std::string &,
                 const std::string)) &
                 adios2::py11::IO::DefineAttribute,
             pybind11::arg("name"), pybind11::arg("ints"), pybind11::arg("variable_name") = "",
             pybind11::arg("separator") = "/", pybind11::return_value_policy::move)

        .def("DefineAttribute",
             (adios2::py11::Attribute(adios2::py11::IO::*)(
                 const std::string &, const std::vector<double> &, const std::string &,
                 const std::string)) &
                 adios2::py11::IO::DefineAttribute,
             pybind11::arg("name"), pybind11::arg("doubles"), pybind11::arg("variable_name") = "",
             pybind11::arg("separator") = "/", pybind11::return_value_policy::move)

        .def("DefineAttribute",
             (adios2::py11::Attribute(adios2::py11::IO::*)(
                 const std::string &, const std::vector<std::complex<double>> &,
                 const std::string &, const std::string)) &
                 adios2::py11::IO::DefineAttribute,
             pybind11::arg("name"), pybind11::arg("complexes"), pybind11::arg("variable_name") = "",
             pybind11::arg("separator") = "/", pybind11::return_value_policy::move)

        .def("DefineAttribute",
             (adios2::py11::Attribute(adios2::py11::IO::*)(
                 const std::string &, const pybind11::object &, const std::string &,
                 const std::string)) &
                 adios2::py11::IO::DefineAttribute,
             pybind11::arg("name"), pybind11::arg("value"), pybind11::arg("variable_name") = "",
             pybind11::arg("separator") = "/", pybind11::return_value_policy::move)

        .def("Open", (adios2::py11::Engine(adios2::py11::IO::*)(const std::string &, const int)) &
                         adios2::py11::IO::Open)
#if ADIOS2_USE_MPI
        .def("Open", (adios2::py11::Engine(adios2::py11::IO::*)(const std::string &, const int,
                                                                adios2::py11::MPI4PY_Comm comm)) &
                         adios2::py11::IO::Open)
#endif
        .def("AvailableAttributes", &adios2::py11::IO::AvailableAttributes,
             pybind11::arg("varname") = "", pybind11::arg("separator") = "/",
             pybind11::return_value_policy::move)

        .def("AvailableVariables", &adios2::py11::IO::AvailableVariables)
        .def("FlushAll", &adios2::py11::IO::FlushAll)
        .def("EngineType", &adios2::py11::IO::EngineType)
        .def("RemoveVariable", &adios2::py11::IO::RemoveVariable)
        .def("RemoveAllVariables", &adios2::py11::IO::RemoveAllVariables)
        .def("RemoveAttribute", &adios2::py11::IO::RemoveAttribute)
        .def("RemoveAllAttributes", &adios2::py11::IO::RemoveAllAttributes);

    pybind11::class_<adios2::py11::Query>(m, "Query")
        .def("__nonzero__",
             [](const adios2::py11::Query &query) {
                 const bool opBool = query ? true : false;
                 return opBool;
             })
        // Python 3
        .def("__bool__",
             [](const adios2::py11::Query &query) {
                 const bool opBool = query ? true : false;
                 return opBool;
             })
        .def(pybind11::init<const std::string &, const adios2::py11::Engine &>(),
             "adios2 query construction, a xml query File and a read engine",
             pybind11::arg("queryFile"), pybind11::arg("reader") = true)

        .def("GetResult", &adios2::py11::Query::GetResult)
        .def("GetBlockIDs", &adios2::py11::Query::GetBlockIDs);

    pybind11::class_<adios2::py11::Operator>(m, "Operator")
        // Python 2
        .def("__nonzero__",
             [](const adios2::py11::Operator &op) {
                 const bool opBool = op ? true : false;
                 return opBool;
             })
        // Python 3
        .def("__bool__",
             [](const adios2::py11::Operator &op) {
                 const bool opBool = op ? true : false;
                 return opBool;
             })
        .def("Type", &adios2::py11::Operator::Type)
        .def("SetParameter", &adios2::py11::Operator::SetParameter)
        .def("Parameters", &adios2::py11::Operator::Parameters);

    pybind11::class_<adios2::py11::Variable>(m, "Variable")
        // Python 2
        .def("__nonzero__",
             [](const adios2::py11::Variable &variable) {
                 const bool opBool = variable ? true : false;
                 return opBool;
             })
        // Python 3
        .def("__bool__",
             [](const adios2::py11::Variable &variable) {
                 const bool opBool = variable ? true : false;
                 return opBool;
             })
        .def("SetShape", &adios2::py11::Variable::SetShape)
        .def("SetBlockSelection", &adios2::py11::Variable::SetBlockSelection)
        .def("SetSelection", &adios2::py11::Variable::SetSelection)
        .def("SetStepSelection", &adios2::py11::Variable::SetStepSelection)
        .def("SelectionSize", &adios2::py11::Variable::SelectionSize)
        .def("Name", &adios2::py11::Variable::Name)
        .def("Type", &adios2::py11::Variable::Type)
        .def("Sizeof", &adios2::py11::Variable::Sizeof)
        .def("ShapeID", &adios2::py11::Variable::ShapeID)
        .def("Shape", &adios2::py11::Variable::Shape,
             pybind11::arg("step") = adios2::EngineCurrentStep)
        .def("Start", &adios2::py11::Variable::Start)
        .def("Count", &adios2::py11::Variable::Count)
        .def("Steps", &adios2::py11::Variable::Steps)
        .def("StepsStart", &adios2::py11::Variable::StepsStart)
        .def("BlockID", &adios2::py11::Variable::BlockID)
        .def("SingleValue", &adios2::py11::Variable::SingleValue)
        .def("AddOperation", (size_t(adios2::py11::Variable::*)(const adios2::py11::Operator,
                                                                const adios2::Params &)) &
                                 adios2::py11::Variable::AddOperation)
        .def("AddOperation",
             (size_t(adios2::py11::Variable::*)(const std::string &, const adios2::Params &)) &
                 adios2::py11::Variable::AddOperation)
        .def("Operations", &adios2::py11::Variable::Operations)
        .def("RemoveOperations", &adios2::py11::Variable::RemoveOperations);

    pybind11::class_<adios2::py11::Attribute>(m, "Attribute")
        // Python 2
        .def("__nonzero__",
             [](const adios2::py11::Attribute &attribute) {
                 const bool opBool = attribute ? true : false;
                 return opBool;
             })
        // Python 3
        .def("__bool__",
             [](const adios2::py11::Attribute &attribute) {
                 const bool opBool = attribute ? true : false;
                 return opBool;
             })
        .def("Name", &adios2::py11::Attribute::Name)
        .def("Type", &adios2::py11::Attribute::Type)
        .def("DataString", &adios2::py11::Attribute::DataString)
        .def("Data", &adios2::py11::Attribute::Data)
        .def("SingleValue", &adios2::py11::Attribute::SingleValue);

    pybind11::class_<adios2::py11::Engine>(m, "Engine")
        // Python 2
        .def("__nonzero__",
             [](const adios2::py11::Engine &engine) {
                 const bool opBool = engine ? true : false;
                 return opBool;
             })
        // Python 3
        .def("__bool__",
             [](const adios2::py11::Engine &engine) {
                 const bool opBool = engine ? true : false;
                 return opBool;
             })
        .def("BeginStep",
             (adios2::StepStatus(adios2::py11::Engine::*)(const adios2::StepMode, const float)) &
                 adios2::py11::Engine::BeginStep,
             pybind11::arg("mode"), pybind11::arg("timeoutSeconds") = -1.f,
             pybind11::return_value_policy::move)

        .def("BeginStep",
             (adios2::StepStatus(adios2::py11::Engine::*)()) & adios2::py11::Engine::BeginStep,
             pybind11::return_value_policy::move)

        .def("Put",
             (void(adios2::py11::Engine::*)(adios2::py11::Variable, const pybind11::array &,
                                            const adios2::Mode launch)) &
                 adios2::py11::Engine::Put,
             pybind11::arg("variable"), pybind11::arg("array"),
             pybind11::arg("launch") = adios2::Mode::Deferred)

        .def("Put", (void(adios2::py11::Engine::*)(adios2::py11::Variable, const std::string &)) &
                        adios2::py11::Engine::Put)

        .def("Put",
             (void(adios2::py11::Engine::*)(adios2::py11::Variable, const std::vector<int64_t> &,
                                            const adios2::Mode launch)) &
                 adios2::py11::Engine::Put,
             pybind11::arg("variable"), pybind11::arg("ints"),
             pybind11::arg("launch") = adios2::Mode::Sync)

        .def("Put",
             (void(adios2::py11::Engine::*)(adios2::py11::Variable, const std::vector<double> &,
                                            const adios2::Mode launch)) &
                 adios2::py11::Engine::Put,
             pybind11::arg("variable"), pybind11::arg("floats"),
             pybind11::arg("launch") = adios2::Mode::Sync)

        .def("Put",
             (void(adios2::py11::Engine::*)(adios2::py11::Variable,
                                            const std::vector<std::complex<double>> &,
                                            const adios2::Mode launch)) &
                 adios2::py11::Engine::Put,
             pybind11::arg("variable"), pybind11::arg("complexes"),
             pybind11::arg("launch") = adios2::Mode::Sync)

        .def("PerformPuts", &adios2::py11::Engine::PerformPuts)

        .def("PerformDataWrite", &adios2::py11::Engine::PerformDataWrite)

        .def("Get",
             (void(adios2::py11::Engine::*)(adios2::py11::Variable, pybind11::array &,
                                            const adios2::Mode launch)) &
                 adios2::py11::Engine::Get,
             pybind11::arg("variable"), pybind11::arg("array"),
             pybind11::arg("launch") = adios2::Mode::Deferred)

        .def("Get",
             (std::string(adios2::py11::Engine::*)(adios2::py11::Variable,
                                                   const adios2::Mode launch)) &
                 adios2::py11::Engine::Get,
             pybind11::arg("variable"), pybind11::arg("launch") = adios2::Mode::Deferred)

        .def("PerformGets", &adios2::py11::Engine::PerformGets)

        .def("EndStep", &adios2::py11::Engine::EndStep)

        .def("BetweenStepPairs", &adios2::py11::Engine::BetweenStepPairs)

        .def("Flush", &adios2::py11::Engine::Flush)

        .def("Close", &adios2::py11::Engine::Close, pybind11::arg("transportIndex") = -1)

        .def("CurrentStep", &adios2::py11::Engine::CurrentStep)

        .def("Name", &adios2::py11::Engine::Name)

        .def("Type", &adios2::py11::Engine::Type)

        .def("Steps", &adios2::py11::Engine::Steps)

        .def("LockWriterDefinitions", &adios2::py11::Engine::LockWriterDefinitions)

        .def("LockReaderSelections", &adios2::py11::Engine::LockReaderSelections)

        .def("BlocksInfo", &adios2::py11::Engine::BlocksInfo);
}
