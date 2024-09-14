## ADIOS2 simulations examples

The _simulations_ examples are meant to demonstrate how to integrate ADIOS2 within your
simulation code to read, write and/or stream your simulation data.

They can be found in the following subdirectories, and they should be explored in the order that they are listed:

1. [korteweg-de-vries](korteweg-de-vries): The _korteweg-de-vries_ example demonstrates how to solve the initial value
   problem for the Korteweg de-Vries equation via the Zabusky and Krustal scheme and integrate ADIOS2 for its IO
   capabilities.
   * Languages: C++, Python
2. [lorenz_ode](lorenz_ode): The _lorenz_ode_ example demonstrates how to solve the Lorenz system of ordinary
   differential equations and integrate ADIOS2 for its IO capabilities.
   * Languages: C++
3. [gray-scott](gray-scott): The _gray-scott_ example demonstrates how to perform simulation with the Gray-Scott
   reaction diffusion model and integrate ADIOS2 for its IO capabilities.
   * Languages: C++
4. [gray-scott-struct](gray-scott-struct): The _gray-scott-struct_ example demonstrates how to perform a simulation with
   the Gray-Scott reaction diffusion model using a struct that defines user defined data types and integrate ADIOS2 for
   its IO capabilities.
   * Languages: C++
5. [gray-scott-kokkos](gray-scott-kokkos): The _gray-scott-kokkos_ example demonstrates how to perform a simulation with
   the Gray-Scott reaction diffusion model using Kokkos and integrate ADIOS2 for its IO capabilities.
   * Languages: C++
6. [GrayScott.jl](GrayScott.jl): The _GrayScott.jl_ example demonstrates how to perform a simulation with the Gray-Scott
   reaction diffusion model and integrate ADIOS2 for its IO capabilities.
   * Languages: Julia
7. [heatTransfer](heatTransfer): The _heatTransfer_ example demonstrates how to solve a 2D Poisson equation for
   temperature in homogeneous media using finite differences and integrate ADIOS2 for its IO capabilities.
   * Languages: C++
