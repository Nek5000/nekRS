**ADIOS2** is the latest implementation of the `Adaptable Input Output System <https://csmd.ornl.gov/software/adios2>`_.
This brand new architecture continues the performance legacy of ADIOS1, and extends its capabilities to address the extreme challenges of scientific data IO.

The `ADIOS2 repo is hosted at GitHub <https://github.com/ornladios/ADIOS2>`_.

The ADIOS2 infrastructure is developed as a multi-institutional collaboration
between

  * `Oak Ridge National Laboratory <https://www.ornl.gov>`_
  * `Kitware Inc. <https://www.kitware.com>`_
  * `Lawrence Berkeley National Laboratory <http://www.lbl.gov>`_
  * `Georgia Institute of Technology <http://www.gatech.edu>`_
  * `Rutgers University <http://www.rutgers.edu>`_

The key aspects ADIOS2 are

#. **Modular architecture:** ADIOS2 takes advantage of the major features
   of C++11. The architecture utilizes a balanced combination of runtime
   polymorphism and template meta-programming to expose intuitive abstractions for a broad range of IO applications.


#. **Community:** By maintaining coding standards, collaborative
   workflows, and understandable documentation, ADIOS2 lowers the barriers to entry for scientists to meaningfully interact with the code.


#. **Sustainability:** Continuous integration and unit testing ensure that ADIOS2 evolves responsibly.
   Bug reports are triaged and fixed in a timely manner and can be reported on `GitHub <https://github.com/ornladios/ADIOS2/issues>`_.


#. **Language Support:** In addition to the native C++, bindings for Python, C, Fortran and Matlab are provided.


#. **Commitment:** ADIOS2 is committed to the HPC community, releasing a new version every 6 months.

*ADIOS2 is funded by the Department of Energy as part of the* `Exascale Computing Project <https://www.exascaleproject.org>`_.

************************
What ADIOS2 is and isn't
************************

**ADIOS2 is:**

- **A Unified High-performance I/O Framework**: using the same abstraction API ADIOS2 can transport and transform groups of self-describing data variables and attributes across different media (file, wide-area-network, in-memory staging, etc.) with performance an ease of use as the main goals.

- **MPI-based**: parallel MPI applications as well as serial codes can use it

- **Streaming-oriented**: ADIOS2 favors codes transferring a group of variables asynchronously wherever possible. Moving one variable at a time, in synchronous fashion, is the special case rather than normal.

- **Step-based**: to resemble actual production of data in "steps" of variable groups, for either streaming or random-access (file) media

- **Free and open-source**: ADIOS2 is permissibly licensed under the OSI-approved Apache 2 license.

- **Extreme scale I/O**: ADIOS2 is being used in supercomputer applications that write and read up to several petabytes in a single simulation run. ADIOS2 is designed to provide scalable I/O on the largest supercomputers in the world.


**ADIOS2 is not**:

- **A file-only I/O library**: Code coupling and in situ analyis is possible through files but special engines are available to achieve the same thing faster through TCP, RDMA and MPI communication. High performance write/read using a file system is a primary goal of ADIOS2 though.

- **MPI-only**

- **A Hierarchical Model**: Data hierarchies can be built on top of the ADIOS2 according to the application, but ADIOS2 sits a layer of abstraction beneath this.

- **A Memory Manager Library**: we don't own or manage the application's memory
