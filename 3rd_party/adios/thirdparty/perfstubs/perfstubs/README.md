# PerfStubs

Copyright (c) 2019-2022 University of Oregon
Distributed under the BSD Software License
(See accompanying file LICENSE.txt)

Profiling API for adding external tool instrumentation support to any project.

This was motivated by the need to quickly add instrumentation to the
[ADIOS2](https://github.com/ornladios/ADIOS2) library without adding a build
dependency, or tying to a specific measurement tool.

The initial prototype implementation was tied to TAU, but evolved to this more
generic version, which was extracted as a separate repository for testing and
demonstration purposes.

For a more extensive readme, see [perfstubs_api/README.md](perfstubs_api/README.md).
