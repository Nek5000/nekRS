# Common Staging Test area

_Note:_ This is a directory designed to hold tests that can be run with a variety of test engines

## Todo Items
- [ ] Make this all vastly simpler while maintaining at least this much functionality

- [ ] Classify tests by what features they depend upon and tag engines with what features they support.  When tests are added, this compatability matrix should tell us what engines it can be run with and have expectations that it will pass.  Automatically run all tests with all engines that they are compatible with.

## Overview

This test suite consists of a core set of tests defined in TestSupp.cmake, where a base test is specified by setting a CMake variable of the form <testname>_CMD.  (See further notes in TestSupp.cmake.)

## Base test suite:
* 1x1 - Simple one reader one writer test, can run MPMD with MPI
* 1x1.NoPreload - 1x1 with Preload explicitly turned off (likely SST specific)
* 1x1.ForcePreload - 1x1 with Preload explicitly turned on (likely SST specific)
* 1x1.Bulk - One to one, with larger data
* 1x1.BulkLockGeometry - One to one, with larger data with geometry explicitly locked on reader and writer.
* 1x1.NoData - A simple test where no data is written, but there are still timesteps
* 2x2.NoData - A simple test where no data is written, but there are still timesteps
* 2x2.HalfNoData - A simple test where one rank writes no data in a timestep, while the other writes all of it, everyone still timesteps
* 2x1 - simple 2 to 1 MPI test
* 2x1.ZeroDataVar - simple 2 to 1 MPI test, with one rank contributing zero data to some variable
* 2x1.NoPreload - 2x1 with Preload explicitly turned off (likely SST specific)
* 2x3.ForcePreload - 2x3 with Preload explicitly turned on (likely SST specific)
* 1x2 - Simple one reader two writer test, can run MPMD with MPI
* 3x5 - Simple 3 reader 5 writer test, can run MPMD with MPI
* 3x5LockGeometry - Simple 3 reader 5 writer test,  with geometry explicitly locked on reader and writer.
* 5x3 - Simple 5 reader 3 writer test, can run MPMD with MPI
* 1x1.Local - simple test with local, not global arrays
* 2x1.Local - simple test with local, not global arrays
* 1x2.Local - simple test with local, not global arrays
* 3x5.Local - simple test with local, not global arrays
* 5x3.Local - simple test with local, not global arrays
* 1x1.LocalVarying - test with local variables, but the dimensions change from rank to rank and timestep to timestep
* 5x3.LocalVarying - test with local variables, but the dimensions change from rank to rank and timestep to timestep
* DelayedReader_3x5 - simple test, but with the reader delayed by 5 seconds to see if the writer waits.  No MPMD.
* FtoC.3x5 - Simple 3 reader 5 writer, but with Fortran writer and C++ reader
* FtoF.3x5 - Simple 3 reader 5 writer, but with Fortran writer and Fortran reader
* 1x1.SharedNothing - Two streams between reader and writer, sharing no ADIOS2 resources
* 1x1.SharedIO - Two streams between reader and writer, both sharing an IO
* 1x1.SharedVar - Two streams between reader and writer, both sharing variables (writing a single variable to both)
* 1x1.SharedNothingSync - Two streams between reader and writer, sharing no ADIOS2 resources using Put Sync
* 1x1.SharedIOSync - Two streams between reader and writer, both sharing an IO using Put Sync
* 1x1.SharedVarSync - Two streams between reader and writer, both sharing variables (writing a single variable to both) using Put Sync
* 2x1.SharedNothing - Two streams between reader and writer, sharing no ADIOS2 resources
* 2x1.SharedIO - Two streams between reader and writer, both sharing an IO
* 2x1.SharedVar - Two streams between reader and writer, both sharing variables (writing a single variable to both)
* 2x1.SharedNothingSync - Two streams between reader and writer, sharing no ADIOS2 resources using Put Sync
* 2x1.SharedIOSync - Two streams between reader and writer, both sharing an IO using Put Sync
* 2x1.SharedVarSync - Two streams between reader and writer, both sharing variables (writing a single variable to both) using Put Sync
* NoReaderNoWait - Tests to see if a writer can run to completion without a reader present (RendezvousReaderCount = 0)
* TimeoutOnOpen - Tests to see if a reader will timeout if no writer ever opens
* 1x1.Modes - Tests to see that we can mix and match Put Sync and Deferred modes and still get good data
* 1x1.Attrs - Tests writing and reading of attributes defined before Open
* FtoC.1x1 - Simple Fortran to C++ test
* CtoF.1x1 - Simple C++ to Fortran test
* FtoF.1x1 - Simple Fortran to Fortran test
* ZFPCompression.1x1 - Simple test using engine-level compression
* ZFPCompression.3x5 - 3x5 test using engine-level compression
* KillReadersSerialized.3x2 - Test to see if readers can be killed mid-stream without disrupting the writer, here at most one reader at a time, killing or starting a reader randomly every 5 seconds.
* KillReadersSerialized.3x6 - Test to see if readers can be killed mid-stream without disrupting the writer, here at most one reader at a time.
* KillReadersSerialized3Max.3x6 - Test to see if readers can be killed mid-stream without disrupting the writer, here at most three readers at a time.
* KillWriter_2x2 - Test to see if a reader gets an appropriate failure when a writer is killed mid-stream
* KillWriterTimeout_2x2 - Test to see if a reader gets an appropriate failure when a writer is killed mid-stream with using BeginStep with timeout
* PreciousTimestep.3x2 - Test to see if new readers always see the "precious" timestep 0, despite arriving late
* PreciousTimestepDiscard.3x2 - Test to see if new readers always see the "precious" timestep 0, despite arriving late, here with a slow reader, so that the writer will discard timesteps due to queueing
* TimeoutReader.1x1 - Testing reader-side BeginStep with timeout
* LatestReader.1x1 - Writer runs faster than reader, reader specifies it wants the LatestTimestep, so we expect to skip some (reader delay after EndStep)
* LatestReaderHold.1x1 - Writer runs faster than reader, reader specifies it wants the LatestTimestep, so we expect to skip some (reader delay between Begin and EndStep)
* DiscardWriter.1x1 - Fast writer, slower reader and a queue policy that should cause timesteps to be discarded.  Are they?

