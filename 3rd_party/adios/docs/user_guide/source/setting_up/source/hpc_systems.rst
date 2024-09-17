.. _HPCBuild:

***********************
Building on HPC Systems
***********************

#. **Modules:** Make sure all "module" dependencies are loaded and that minimum requirements are satisfied.
   Load the latest CMake module as many HPC systems default to an outdated version.
   Build with a C++11-compliant compiler, such as gcc >= 4.8.1, Intel >= 15, and PGI >= 15.

#. **Static/Dynamic build:** On Cray systems such as `Titan <https://www.olcf.ornl.gov/kb_articles/compiling-and-node-types/>`_,
   the default behavior is static linkage, thus CMake builds ADIOS2 creates the static library ``libadios2.a`` by default.
   Read the system documentation to enable dynamic compilation, usually by setting an environment variable such as ``CRAYPE_LINK_TYPE=dynamic``.
   Click `here <https://github.com/ornladios/ADIOS2/tree/master/scripts/runconf/runconf_olcf.sh>`_ for a fully configurable script example on OLCF systems.

#. **Big Endian and 32-bit systems:** ADIOS2 hasn't been tested on big
   endian and generally will not build on 32-bit systems. Please be
   aware before attempting to run. 

#. **PGI compilers and C++11 support:** Version 15 of the PGI compiler is C++11 compliant.
   However it relies on the C++ standard library headers supplied by the system version of GCC, which may or may support all the C++11 features used in ADIOS2.
   On many systems (Titan at OLCF, for example) even though the PGI compiler supports C++11, the configured GCC and its headers do not (4.3.x on Cray Linux Environment, and v5 systems like Titan).
   To configure the PGI compiler to use a newer GCC, you must create a configuration file in your home directory that overrides the PGI compiler's default configuration.
   On Titan, the following steps will re-configure the PGI compiler to use GCC 6.3.0 provided by a module:

.. code-block:: bash

  $ module load gcc/6.3.0
  $ makelocalrc $(dirname $(which pgc++)) -gcc $(which gcc) -gpp $(which g++) -g77 $(which gfortran) -o -net 1>${HOME}/.mypgirc 2>/dev/null



#. **Enabling RDMA for SST data transfers:** The SST engine in ADIOS2
   is capable of using RDMA networks for transfering data between
   writer and reader cohorts, and generally this is the most
   performant data transport.  However, SST depends upon libfabric to
   provide a generic interface to the underlying RDMA capabilities of
   the network, and properly configuring libfabric can be a difficult
   and error-prone task.  HPC computing resources tend to be one-off
   custom resources with their own idiosyncracies, so this
   documentation cannot offer a definitive guide for every situation,
   but we can provide some general guidance and some recommendations
   for specific machines.  If you are unable to configure ADIOS2 and
   libfabric to use RDMA, the best way to get help is to open an issue
   on the ADIOS2 github repository. 

Pre-build concerns of note:

   	* on some HPC resources, libfabric is available as a loadable
	  module.  That should not be taken as an indication that that
	  build of libfabric will work with SST, or even that it is
	  compatible with the system upon which you find it.  Your
	  mileage may vary and you may have to build libfabric
	  manually. 
	* libfabric itself depends upon other libraries like
	  libibverbs and librdmacm.  If you build libfabric with a
	  package manager like spack, spack may build custom versions
	  of those libraries as well, which may conflict with the
	  system versions of those libraries. 
	* MPI on your HPC system may use libfabric itself, and linking
	  your application with a different version of libfabric (or
	  its dependent libraries) may result failure, possibly
	  including opaque error messages from MPI. 
	* libfabric is structured in such a way that even if it is
	  found during configuration, ADIOS *cannot* determine at
	  compile time what providers will be present at run-time, or
	  what their capabilities are.  Therefore even a build that
	  seems to successfully include libfabric and RDMA may be
	  rejected at runtime as unable to support SST data transfer. 

Configuration:
	ADIOS2 uses the CMake find_package() functionality to locate
	libfabric.  CMake will automatically search system libraries,
	but if you need to specify a libfabric location other than in
	a default system location you can add a
	"-DLIBFABRIC_ROOT=<directory>" argument to direct CMake to
	libfabric's location.   If CMake finds libfabric, you should
	see the line "RDMA Transport for Staging: Available" near the
	end of the CMake output. This makes the RDMA DataTransport the
	default for SST data movement.  (More information about SST
	engine parameters like `DataTransport` appears in the SST
	engine description.)  If instead you see "RDMA Transport for
	Staging: Unconfigured", RDMA will not be available to SST.

Run-time:
	Generally, if RDMA is configured and the libfabric provider
	has the capabilities that SST needs for RDMA data transfer,
	SST will use RDMA without external acknowledgement.  However,
	if RDMA is configured, but the libfabric provider doesn't have
	the capabilities that SST needs, ADIOS will output an error :
	'Warning:  Preferred DataPlane "RDMA" not found.'
	If you see this warning in a situation where you expect RDMA
	to be used, enabling verbose debugging output from SST may
	provide more information.  The SstVerbose environment
	variable can have values from 1 to 5, with 1 being minimal
	debugging info (such as confirming which DataTransport is
	being used), and 5 being the most detailed debugging
	information from all ranks.
