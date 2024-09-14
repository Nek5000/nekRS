/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Solves the initial value problem for the Korteweg de-Vries equation via the
 * Zabusky and Krustal scheme.
 *
 * See: https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.15.240 Zabusky, Norman
 * J., and Martin D. Kruskal. "Interaction of solitons in a collisionless plasma
 * and the recurrence of initial states." Physical review letters 15.6 (1965): 240.
 */
#include <adios2.h>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

void display_progress(double progress)
{
    int barWidth = 70;

    std::cout << "\033[0;32m[";
    int pos = static_cast<int>(barWidth * progress);
    for (int i = 0; i < barWidth; ++i)
    {
        if (i < pos)
            std::cout << "=";
        else if (i == pos)
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << "%\033[0m\r";
    std::cout.flush();
}

// This momentum is a conserved quantity of the numerical method.
// Use it to sanity check the solution.
template <typename Real>
Real momentum(std::vector<Real> const &u)
{
    Real p = 0;
    for (size_t i = 0; i < u.size(); ++i)
    {
        p += u[i];
    }
    return p;
}

template <typename Real>
void KdV(int64_t N, Real dt, Real t_max, Real delta = 0.022)
{
    using std::cos;
    const Real pi = 4 * std::atan(Real(1));

    if (N <= 0)
    {
        throw std::domain_error("N > 0 is required");
    }
    if (dt > 1)
    {
        throw std::domain_error("time step is too big");
    }
    if (dt <= 0)
    {
        throw std::domain_error("dt > 0 is required");
    }

    Real dx = static_cast<Real>(2) / (static_cast<Real>(N));

    int64_t M = static_cast<int64_t>(std::ceil(t_max / dt));
    std::cout << "Solving the initial value problem for the KdV equation ∂tu + "
                 "u∂ₓu + δ²∂ₓ³u = 0 using δ = "
              << delta << ".\n";
    std::cout << "Initial conditions: u(x,0) = cos(πx) for xϵ[0,2].\n";
    std::cout << "Using ∆x = " << dx << ", ∆t = " << dt << " and t_max = " << t_max << "\n";
    adios2::ADIOS adios;
    adios2::IO io = adios.DeclareIO("myio");
    auto u_variable =
        io.DefineVariable<Real>("u", {size_t(N)}, {0}, {size_t(N)}, adios2::ConstantDims);
    io.DefineAttribute<Real>("x0", 0);
    io.DefineAttribute<Real>("dx", dx);
    io.DefineAttribute<std::string>("interpretation", "Equispaced");
    adios2::Engine adios_engine = io.Open("korteweg_de_vries.bp", adios2::Mode::Write);

    std::vector<Real> u0(N);
    for (int64_t i = 0; i < N; ++i)
    {
        u0[i] = cos(pi * static_cast<double>(i) * dx);
    }

    adios_engine.BeginStep();
    adios_engine.Put(u_variable, u0.data());
    adios_engine.EndStep();
    // Real original_momentum = momentum(u0);
    std::vector<Real> u1(N);
    for (int64_t i = 0; i < N; ++i)
    {
        Real cdt = cos(pi * static_cast<double>(i) * dx) * dt;
        u1[i] = cos(pi * (static_cast<double>(i) * dx - cdt));
    }

    Real k1 = dt / (3 * dx);
    Real k2 = delta * delta * dt / (dx * dx * dx);
    Real t1, t2;
    std::vector<Real> u2(N);
    for (int64_t j = 1; j < M - 1; ++j)
    {
        t1 = (u1[1] + u1[0] + u1[N - 1]) * (u1[1] - u1[N - 1]);
        t2 = u1[2] - 2 * u1[1] + 2 * u1[N - 1] - u1[N - 2];
        u2[0] = u0[0] - k1 * t1 - k2 * t2;
        t1 = (u1[2] + u1[1] + u1[0]) * (u1[2] - u1[0]);
        t2 = u1[3] - 2 * u1[2] + 2 * u1[0] - u1[N - 1];
        u2[1] = u0[1] - k1 * t1 - k2 * t2;
        for (int64_t i = 2; i < N - 2; ++i)
        {
            t1 = (u1[i + 1] + u1[i] + u1[i - 1]) * (u1[i + 1] - u1[i - 1]);
            t2 = u1[i + 2] - 2 * u1[i + 1] + 2 * u1[i - 1] - u1[i - 2];
            u2[i] = u0[i] - k1 * t1 - k2 * t2;
        }
        u2[N - 2] = u0[N - 2] - k1 * (u1[N - 1] + u1[N - 2] + u1[N - 3]) * (u1[N - 1] - u1[N - 3]) -
                    k2 * (u1[0] - 2 * u1[N - 1] + 2 * u1[N - 3] - u1[N - 4]);
        u2[N - 1] = u0[N - 1] - k1 * (u1[0] + u1[N - 1] + u1[N - 2]) * (u1[0] - u1[N - 2]) -
                    k2 * (u1[1] - 2 * u1[0] + 2 * u1[N - 2] - u1[N - 3]);

        Real p = momentum(u2);
#ifndef __INTEL_LLVM_COMPILER
        if (std::isnan(p))
        {
            std::cerr << "\nSolution diverged at t = " << (static_cast<double>(j + 1)) * dt << "\n";
            std::cerr << "Momentum = " << p << "\n";
            return;
        }
#endif
        if (std::abs(p) > std::sqrt(std::numeric_limits<Real>::epsilon()))
        {
            std::cerr << "\nSolution diverged at t = " << (static_cast<double>(j + 1)) * dt << "\n";
            std::cerr << "Momentum = " << p << "\n";
        }
        // The stability condition is very severe this iteration: We have to
        // take way more time steps for stability than we need for accuracy.
        // Hence we need skip a ton of steps so we don't spend all our time
        // writing data:
        int64_t skip_steps = 40000;
        if ((j + 1) % skip_steps == 0)
        {
            display_progress(static_cast<double>(j + 1) / (static_cast<double>(M - 1)));
            adios_engine.BeginStep();
            adios_engine.Put(u_variable, u2.data());
            adios_engine.EndStep();
        }
        u0 = u1;
        u1 = u2;
    }
    adios_engine.Close();
}

int main(int argc, char **argv)
{
    int64_t N = 256;
    double t_max = 5;
    double delta = 0.022;
    if (argc > 1)
    {
        std::string dx_str = argv[1];
        if (dx_str == "-h" || dx_str == "--help")
        {
            std::cout << "Usage: ./adios2_simulations_kortewegDeVries"
                         "N t_max δ, where N is number of spacial "
                         "gridpoints (∆t chosen from ∆x via the stability "
                         "condition), t_max is max simulation time, and δ is "
                         "an interaction parameter;"
                         "e.g., ./adios2_simulations_kortewegDeVries 512 10 0.022\n";
            return 0;
        }
        N = std::stoi(dx_str);
    }
    if (argc > 2)
    {
        t_max = std::stod(argv[2]);
    }
    if (argc > 3)
    {
        delta = std::stod(argv[3]);
    }

    double dx = 1.0 / static_cast<double>(N);
    double dt = 27 * dx * dx * dx / 4;
    try
    {
        KdV<double>(N, dt, t_max, delta);
    }
    catch (std::exception const &e)
    {
        std::cerr << "Caught exception from KdV call: " << e.what() << "\n";
        return 1;
    }
}
