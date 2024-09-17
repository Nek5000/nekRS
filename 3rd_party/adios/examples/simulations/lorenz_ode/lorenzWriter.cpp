/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */

#include "lorenz.hpp"
#include <adios2.h>
#include <array>
#include <cassert>
#include <iostream>
#include <vector>

void solve_lorenz_ivp()
{
    using Real = double;
    Real sigma = 10;
    Real beta = Real(8) / Real(3);
    Real rho = 28;
    Real tmax = 10;
    Real absolute_error = 1e-5;

    adios2::ADIOS adios;
    adios2::IO io = adios.DeclareIO("myio");
    io.DefineAttribute<Real>("σ", sigma);
    io.DefineAttribute<Real>("β", beta);
    io.DefineAttribute<Real>("ρ", rho);
    io.DefineAttribute<Real>("‖û-u‖", absolute_error);
    io.DefineAttribute<std::string>("interpretation",
                                    "7D array of structs {tᵢ, xᵢ, yᵢ, zᵢ, ẋᵢ, ẏᵢ, żᵢ}");
    adios2::Engine adios_engine = io.Open("lorenz.bp", adios2::Mode::Write);
    // Check that std::vector<std::array<Real, 7>> has the proper memory layout.
    static_assert(sizeof(std::array<Real, 7>) == 7 * sizeof(Real),
                  "The std::array on your system does not have the proper "
                  "layout to be correctly serialized in ADIOS2.");

    size_t paths = 2;
    for (size_t i = 0; i < paths; ++i)
    {
        for (size_t j = 0; j < paths; ++j)
        {
            for (size_t k = 0; k < paths; ++k)
            {
                const std::array<Real, 3> initial_conditions{Real(i), Real(j), Real(k)};
                auto solution =
                    lorenz<double>(sigma, beta, rho, initial_conditions, tmax, absolute_error);
                auto const &skeleton = solution.states();
                // This is just being *ultra* cautious; this double checks that
                // the memory layout is contiguous.
                const double *const p = skeleton[0].data();
                for (size_t l = 0; l < skeleton.size(); ++l)
                {
                    for (size_t m = 0; m < 7; ++m)
                    {
                        assert(skeleton[l][m] == p[7 * l + m]);
                    }
                }
                const std::string variable =
                    "u" + std::to_string(i) + std::to_string(j) + std::to_string(k);
                auto state_variable =
                    io.DefineVariable<Real>(variable, {7 * skeleton.size()}, {0},
                                            {7 * skeleton.size()}, adios2::ConstantDims);
                adios_engine.Put(state_variable, p, adios2::Mode::Sync);
            }
        }
    }
    adios_engine.Close();
}

int main()
{
    try
    {
        solve_lorenz_ivp();
    }
    catch (std::exception const &e)
    {
        std::cerr << "Caught an exception solving Lorenz ODE: " << e.what() << "\n";
    }

    try
    {
        test_lorenz<double>();
    }
    catch (std::exception const &e)
    {
        std::cout << "Caught and exception from the Lorenz unit tests: " << e.what() << "\n";
    }
}
