[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Documentation](https://readthedocs.org/projects/adios2/badge/?version=latest)](https://adios2.readthedocs.io/en/latest/?badge=latest)

[![GitHub (pre-)release](https://img.shields.io/github/release/ornladios/adios2/all.svg)]()
[![Spack Version](https://img.shields.io/spack/v/adios2.svg)](https://spack.readthedocs.io/en/latest/package_list.html#adios2)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/adios2)](https://anaconda.org/conda-forge/adios2)
[![PyPI version](https://badge.fury.io/py/adios2.svg)](https://badge.fury.io/py/adios2)

[![Circle CI](https://circleci.com/gh/ornladios/ADIOS2.svg?style=shield)](https://circleci.com/gh/ornladios/ADIOS2)

[![Coverity Scan Build Status](https://scan.coverity.com/projects/11116/badge.svg)](https://scan.coverity.com/projects/ornladios-adios2)


# ADIOS2 : The Adaptable Input Output System version 2

This is ADIOS2: The Adaptable Input/Output (I/O) System.

ADIOS2 is developed as part of the United States Department of Energy's Exascale Computing Project.
It is a framework for scientific data I/O to publish and subscribe to data when and where required.

ADIOS2 transports data as groups of self-describing variables and attributes across different media types (such as files, wide-area-networks, and remote direct memory access) using a common application programming interface for all transport modes.
ADIOS2 can be used on supercomputers, cloud systems, and personal computers.

ADIOS2 focuses on:

1. **Performance** I/O scalability in high performance computing (HPC) applications.
2. **Adaptability** unified interfaces to allow for several modes of transport (files, memory-to-memory)
3. **Ease of Use** two-level application programming interface (APIs)
* Full APIs for HPC applications: C++11, Fortran 90, C 99, Python 2 and 3
* Simplified High-Level APIs for data analysis: Python 2 and 3, C++11, Matlab

In addition, ADIOS2 APIs are based on:

* **MPI** Although ADIOS2 is MPI-based, it can also be used in non-MPI serial code.

* **Data Groups** ADIOS2 favors a deferred/prefetch/grouped variables transport mode by default to maximize data-per-request ratios.
Sync mode, one variable at a time, is treated as the special case.

* **Data Steps** ADIOS2 follows the actual production/consumption of data using an I/O “steps” abstraction removing the need to manage extra indexing information.

* **Data Engines** ADIOS2 Engine abstraction allows for reusing the APIs for different transport modes removing the need for drastic code changes.

## Documentation

Documentation is hosted at [readthedocs](https://adios2.readthedocs.io).

## Citing

If you find ADIOS2 useful, please cite our [SoftwareX paper](https://doi.org/10.1016/j.softx.2020.100561), which also gives a high-level overview to the motivation and goals of ADIOS; complementing the documentation.

## Getting ADIOS2

* From source: [Install ADIOS2 documentation](https://adios2.readthedocs.io/en/latest/setting_up/setting_up.html#).
For a `cmake` configuration example see [scripts/runconf/runconf.sh](https://github.com/ornladios/ADIOS2/blob/master/scripts/runconf/runconf.sh)


* Conda packages:
* [https://anaconda.org/conda-forge/adios2](https://anaconda.org/conda-forge/adios2)


* Spack: [adios2 package](https://spack.readthedocs.io/en/latest/package_list.html#adios2)


* Docker images: under [scripts/docker](https://github.com/ornladios/ADIOS2/tree/master/scripts/docker)


Once ADIOS2 is installed refer to:

* [Linking ADIOS2](https://adios2.readthedocs.io/en/latest/setting_up/setting_up.html#linking-adios-2)


## Releases

* Latest release: [v2.10.1](https://github.com/ornladios/ADIOS2/releases/tag/v2.10.1)

* Previous releases: [https://github.com/ornladios/ADIOS2/releases](https://github.com/ornladios/ADIOS2/releases)

## Community

ADIOS2 is an open source project: Questions, discussion, and contributions are welcome. Join us at:

- Mailing list: adios-ecp@kitware.com
- Github Discussions: https://github.com/ornladios/ADIOS2/discussions

## Reporting Bugs

If you find a bug, please open an [issue on ADIOS2 github repository](https://github.com/ornladios/ADIOS2/issues)

## Contributing

See the [Contributor's Guide to ADIOS 2](Contributing.md) for instructions on how to contribute.

## License
ADIOS2 is licensed under the Apache License v2.0.
See the accompanying [Copyright.txt](Copyright.txt) for more details.

## Directory layout

* bindings - public application programming interface, API, language bindings (C++11, C, Fortran, Python and Matlab)

* cmake - Project specific CMake modules

* examples - Simple set of examples in different languages

* scripts - Project maintenance and development scripts

* source - Internal source code for private components
* adios2 - source directory for the ADIOS2 library to be installed under install-dir/lib/libadios2.
* utils  - source directory for the binary utilities, to be installed under install-dir/bin

* testing - Tests using [gtest](https://github.com/google/googletest)
