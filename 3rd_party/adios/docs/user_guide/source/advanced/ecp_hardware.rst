######################
ADIOS2 in ECP hardware
######################

ADIOS2 is widely used in ECP (Exascale Computing Project) HPC (high performance
computing) systems, some particular ADIOS2 features needs from specifics
workarounds to run successfully.

OLCF CRUSHER
============

SST MPI Data Transport
----------------------

MPI Data Transport relies on client-server features of MPI which are currently
supported in Cray-MPI implementations with some caveats. Here are some of the
observed issues and what its workaround (if any) are:

**MPI_Finalize** will block the system process in the "Writer/Producer" ADIOS2
instance. The reason is that the Producer ADIOS instance internally calls
`MPI_Open_port` which somehow even after calling `MPI_Close_port` `MPI_Finalize`
still consider its port to be in used, hence blocking the process. The
workaround is to use a `MPI_Barrier(MPI_COMM_WORLD)` instead of `MPI_Finalize()`
call.

**srun does not understand mpmd instructions** Simply disable them with the flag
`-DADIOS2_RUN_MPI_MPMD_TESTS=OFF`

**Tests timeout**  Since we launch every tests with srun the scheduling times
can exceed the test default timeout. Use a large timeout (5mins) for running
your tests.

Examples of launching ADIOS2 SST unit tests using MPI DP:

.. code-block:: bash

  # We omit some of the srun (SLURM) arguments which are specific of the project
  # you are working on. Note that you could avoid calling srun directly by
  # setting the CMAKE variable `MPIEXEC_EXECUTABLE`.

  # Launch simple writer test instance
  srun {PROJFLAGS } -N 1 /gpfs/alpine/proj-shared/csc331/vbolea/ADIOS2-build/bin/TestCommonWrite SST mpi_dp_test CPCommPattern=Min,MarshalMethod=BP5

  # On another terminal launch multiple instances of the Reader test
  srun {PROJFLAGS} -N 2 /gpfs/alpine/proj-shared/csc331/vbolea/ADIOS2-build/bin/TestCommonRead SST mpi_dp_test

Alternatively, you can configure your CMake build to use srun directly:

.. code-block:: bash

  cmake . -DMPIEXEC_EXECUTABLE:FILEPATH="/usr/bin/srun" \
       -DMPIEXEC_EXTRA_FLAGS:STRING="-A{YourProject} -pbatch -t10" \
       -DMPIEXEC_NUMPROC_FLAG:STRING="-N" \
       -DMPIEXEC_MAX_NUMPROCS:STRING="-8" \
       -DADIOS2_RUN_MPI_MPMD_TESTS=OFF

  cmake --build .
  ctest

  # monitor your jobs
  watch -n1 squeue -l -u $USER
