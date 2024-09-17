/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * InlineIO.cpp
 *
 *  Created on: May 2020
 *      Author: Caitlin Ross
 */

#include "InlineIO.h"

InlineIO::InlineIO(const Settings &s, MPI_Comm comm)
{
    std::string writerName = "writer";
    std::string readerName = "reader";

    ad = adios2::ADIOS(s.configfile, comm);
    inlineIO = ad.DeclareIO("inlineIO");

    auto params = inlineIO.Parameters();
    if (params.find("writerID") != params.end())
    {
        writerName = params["writerID"];
    }
    if (params.find("readerID") != params.end())
    {
        readerName = params["readerID"];
    }

    if (!inlineIO.InConfigFile())
    {
        inlineIO.SetEngine("Inline");
        inlineIO.SetParameters({{"writerID", writerName}, {"readerID", readerName}});
    }

    // define T as 2D global array
    varT =
        inlineIO.DefineVariable<double>("T",
                                        // Global dimensions
                                        {s.gndx, s.gndy},
                                        // starting offset of the local array in the global space
                                        {s.offsx, s.offsy},
                                        // local size, could be defined later using SetSelection()
                                        {s.ndx, s.ndy});

    inlineWriter = inlineIO.Open(writerName, adios2::Mode::Write, comm);
    inlineReader = inlineIO.Open(readerName, adios2::Mode::Read, comm);

    // Promise that we are not going to change the variable sizes nor add new
    // variables
    inlineWriter.LockWriterDefinitions();
    inlineReader.LockReaderSelections();
}

InlineIO::~InlineIO()
{
    inlineWriter.Close();
    inlineReader.Close();
}

void InlineIO::write(const HeatTransfer &ht)
{
    inlineWriter.BeginStep();
    v = ht.data_noghost();
    inlineWriter.Put<double>(varT, v.data());
    inlineWriter.EndStep();
}

const double *InlineIO::read(bool firstStep)
{
    adios2::StepStatus status = inlineReader.BeginStep(adios2::StepMode::Read);
    if (status != adios2::StepStatus::OK)
    {
        return nullptr;
    }

    auto blocksInfo = inlineReader.BlocksInfo(varT, inlineReader.CurrentStep());
    // in this example we're only expecting one block
    if (blocksInfo.size() != 1)
    {
        throw std::runtime_error("InlineIO::read found incorrect number of blocks");
    }

    auto &info = blocksInfo[0];
    varT.SetBlockSelection(info.BlockID);
    inlineReader.Get(varT, info);
    inlineReader.EndStep();

    // now we can get the pointer with info.Data()
    return info.Data();
}
