/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IO_h5mixer.cpp
 *
 *  Created on: Feb 2017
 *      Author: Norbert Podhorszki
 */

#include "IO.h"

#include <string>

#include <adios2.h>

#define str_helper(X) #X
#define str(X) str_helper(X)
// #ifndef DEFAULT_CONFIG
// #define DEFAULT_CONFIG config.xml
// #endif
#define DEFAULT_CONFIG mix.xml
#define DEFAULT_CONFIG_STR str(DEFAULT_CONFIG)

adios2::ADIOS *ad = nullptr;
// std::shared_ptr<adios2::Engine> h5mixerWriter;
adios2::Engine h5mixerWriter;
adios2::Variable<double> varT;
adios2::Variable<unsigned int> varGndx;

IO::IO(const Settings &s, MPI_Comm comm)
{
    m_outputfilename = MakeFilename(s.outputfile, ".h5");

    /*ad = new adios2::ADIOS(std::string(DEFAULT_CONFIG_STR), comm);
     */
    ad = new adios2::ADIOS(comm);

    // Define method for engine creation

    adios2::IO h5io = ad->DeclareIO("writer");
    if (!h5io.InConfigFile())
    {
        // if not defined by user, we can change the default settings
        // BPFile is the default engine

        // Allow an extra thread for data processing
        // ISO-POSIX file is the default transport
        // Passing parameters to the transport
        h5io.SetEngine("HDFMixer");
    }

    varGndx = h5io.DefineVariable<unsigned int>("gndx");
    h5io.DefineVariable<unsigned int>("gndy");

    // define T as 2D global array
    varT = h5io.DefineVariable<double>("T",
                                       // Global dimensions
                                       {s.gndx, s.gndy},
                                       // starting offset of the local array in the global space
                                       {s.offsx, s.offsy},
                                       // local size, could be defined later using SetSelection()
                                       {s.ndx, s.ndy});

    // add transform to variable
    // adios2::Transform tr = adios2::transform::BZIP2( );
    // varT.AddTransform( tr, "" );
    // varT.AddTransform( tr,"accuracy=0.001" );  // for ZFP

    h5mixerWriter = h5io.Open(m_outputfilename, adios2::Mode::Write, comm);

    if (!h5mixerWriter)
    {
        throw std::ios_base::failure("ERROR: failed to open ADIOS h5mixerWriter\n");
    }
}

IO::~IO()
{
    h5mixerWriter.Close();
    delete ad;
}

void IO::write(int step, const HeatTransfer &ht, const Settings &s, MPI_Comm comm)
{

    h5mixerWriter.BeginStep();
    /* This selection is redundant and not required, since we defined
     * the selection already in DefineVariable(). It is here just as an example.
     */
    // Make a selection to describe the local dimensions of the variable we
    // write and its offsets in the global spaces. This could have been done in
    // adios.DefineVariable()
    varT.SetSelection(adios2::Box<adios2::Dims>({s.offsx, s.offsy}, {s.ndx, s.ndy}));

    /* Select the area that we want to write from the data pointer we pass to
         the
         writer.
         Think HDF5 memspace, just not hyperslabs, only a bounding box
       selection. Engine will copy this bounding box from the data pointer into
       the output buffer. Size of the bounding box should match the "space"
       selection which was given above. Default memspace is always the full
       selection.
    */
    // varT.SetMemorySelection(adios2::Box<adios2::Dims>({1, 1}, {s.ndx,
    // s.ndy}));

    h5mixerWriter.Put<unsigned int>(varGndx, s.gndx);
    h5mixerWriter.Put<unsigned int>("gndy", s.gndy);
    h5mixerWriter.Put<double>(varT, ht.data_noghost().data());

    h5mixerWriter.EndStep();

#ifdef NEVER
#if 1

    /* This selection is redundant and not required, since we defined
     * the selection already in DefineVariable(). It is here just as an example.
     */
    // Make a selection to describe the local dimensions of the variable we
    // write and its offsets in the global spaces. This could have been done in
    // adios.DefineVariable()
    adios2::SelectionBoundingBox sel({s.offsx, s.offsy}, {s.ndx, s.ndy});
    varT->SetSelection(sel);

    /* Select the area that we want to write from the data pointer we pass to
       the
       writer.
       Think HDF5 memspace, just not hyperslabs, only a bounding box selection.
       Engine will copy this bounding box from the data pointer into the output
       buffer.
       Size of the bounding box should match the "space" selection which was
       given
       above.
       Default memspace is always the full selection.
    */
    adios2::SelectionBoundingBox memspace = adios2::SelectionBoundingBox({1, 1}, {s.ndx, s.ndy});
    varT->SetMemorySelection(memspace);

    h5mixerWriter->Write<unsigned int>(*varGndx, s.gndx);
    h5mixerWriter->Write<unsigned int>("gndy", s.gndy);
    h5mixerWriter->Write<double>(*varT, ht.data_noghost().data());
    h5mixerWriter->Advance();

#else
    h5mixerWriter->Write<double>(*varT, ht.data_noghost().data());
    h5mixerWriter->Advance();
#endif
#endif
}
