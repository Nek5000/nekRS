/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */

#ifndef LORENZ_HPP
#define LORENZ_HPP
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

// Solves the Lorenz system using Taylor methods.
// See: https://en.wikipedia.org/wiki/Lorenz_system
// and Corliss, "A Graduate Introduction to Numerical Methods"
// (https://doi.org/10.1007/978-1-4614-8453-0) for an introduction to Taylor
// methods for ODEs.

template <typename Real>
class lorenz
{
public:
    lorenz(const Real sigma, const Real beta, const Real rho,
           std::array<Real, 3> const &initial_conditions, const Real tmax,
           const Real absolute_error_goal)
    {
        using std::sqrt;
        using std::abs;
        using std::cbrt;
        if (tmax <= 0)
        {
            throw std::domain_error("tmax > 0 is required");
        }
        if (absolute_error_goal <= std::numeric_limits<Real>::epsilon())
        {
            throw std::domain_error("Abolute error goal > eps is required");
        }
        Real t = 0;
        Real x = initial_conditions[0];
        Real y = initial_conditions[1];
        Real z = initial_conditions[2];
        Real dotx = sigma * (y - x);
        Real doty = x * (rho - z) - y;
        Real dotz = x * y - beta * z;
        Real ddotx = sigma * (doty - dotx);
        Real ddoty = dotx * (rho - z) - x * dotz - doty;
        Real ddotz = dotx * y + x * doty - beta * dotz;
        states_.emplace_back(std::array<Real, 7>{t, x, y, z, dotx, doty, dotz});
        while (t < tmax)
        {
            // ∆t must satisfy three constraints:
            // ∆t^2 < 2µ/|ddot(x)|, ∆t^2 < 2µ/|ddot(y)|, ∆t^2 < 2µ/|ddot(z)|.
            // where µ = absolute_error_goal.
            const Real m = std::min({abs(ddotx), abs(ddoty), abs(ddotz)});
            Real dt;
            // If all second derivaties are zero, we're actually *more*
            // accurate:
            if (m == 0)
            {
                dt = cbrt(6 * absolute_error_goal);
            }
            else
            {
                dt = sqrt(2 * absolute_error_goal / m);
            }

            // Taylor series:
            t += dt;
            x += dt * dotx + dt * dt * ddotx / 2;
            y += dt * doty + dt * dt * ddoty / 2;
            z += dt * dotz + dt * dt * ddotz / 2;

            // Now compute the derivatives at the new location:
            dotx = sigma * (y - x);
            doty = x * (rho - z) - y;
            dotz = x * y - beta * z;
            ddotx = sigma * (doty - dotx);
            ddoty = dotx * (rho - z) - x * dotz - doty;
            ddotz = dotx * y + x * doty - beta * dotz;

            // And store the state:
            states_.emplace_back(std::array<Real, 7>{t, x, y, z, dotx, doty, dotz});
        }
    }

    // Load data:
    lorenz(std::vector<std::array<Real, 7>> &&state) : states_{std::move(state)}
    {
        // Simple validation: The times increase:
        Real t = states_[0][0];
        // t0 = 0: This is just an incidental feature of the solution,
        // obviously we could change this so that t0 could be arbitrary.
        // But for now, t0 is not arbitrary, so let's use this to validate the
        // deserialization:
        if (t != 0)
        {
            throw std::logic_error("t0 != 0");
        }
        for (size_t i = 1; i < states_.size(); ++i)
        {
            Real ti = states_[i][0];
            if (ti <= t)
            {
                throw std::logic_error("Deserialization is incorrect: Times are not sorted in "
                                       "increasing order t_0 < t_1 < ...");
            }
            t = ti;
        }
    }

    std::array<Real, 3> operator()(Real t) const
    {
        if (t > tmax() || t < tmin())
        {
            throw std::domain_error("t is not in domain of interpolation.");
        }
        if (t == tmax())
        {
            auto const &state = states_.back();
            return {state[1], state[2], state[3]};
        }
        auto comparator = [&](const Real t, std::array<Real, 7> const &state) {
            return t < state[0];
        };
        auto it = std::upper_bound(states_.begin(), states_.end(), t, comparator);
        auto i = std::distance(states_.begin(), it) - 1;
        auto const &s0 = states_[i];
        auto const &s1 = states_[i + 1];
        Real t0 = s0[0];
        Real x0 = s0[1];
        Real y0 = s0[2];
        Real z0 = s0[3];
        Real dx0dt = s0[4];
        Real dy0dt = s0[5];
        Real dz0dt = s0[6];

        Real t1 = s1[0];
        Real x1 = s1[1];
        Real y1 = s1[2];
        Real z1 = s1[3];
        Real dxdt1 = s1[4];
        Real dydt1 = s1[5];
        Real dzdt1 = s1[6];
        // Map t into [0,1]:
        Real dt = s1[0] - s0[0];
        Real s = (t - t0) / dt;
        Real x = (1 - s) * (1 - s) * (x0 * (1 + 2 * s) + dx0dt * (t - t0)) +
                 s * s * (x1 * (3 - 2 * s) + dt * dxdt1 * (s - 1));
        Real y = (1 - s) * (1 - s) * (y0 * (1 + 2 * s) + dy0dt * (t - t0)) +
                 s * s * (y1 * (3 - 2 * s) + dt * dydt1 * (s - 1));
        Real z = (1 - s) * (1 - s) * (z0 * (1 + 2 * s) + dz0dt * (t - t0)) +
                 s * s * (z1 * (3 - 2 * s) + dt * dzdt1 * (s - 1));
        return {x, y, z};
    }

    const std::vector<std::array<Real, 7>> &states() const { return states_; }

    Real tmax() const { return states_.back()[0]; }

    Real tmin() const { return states_.front()[0]; }

    friend std::ostream &operator<<(std::ostream &out, lorenz const &l)
    {
        for (auto &state : l.states_)
        {
            Real t = state[0];
            out << "u(" << t << ") = {" << state[1] << ", " << state[2] << ", " << state[3]
                << "}\n";
        }
        return out;
    }

private:
    std::vector<std::array<Real, 7>> states_;
};

template <typename Real>
void test_lorenz()
{
    using std::abs;
    // Kinda ridiculous to run tests this way,
    // but I don't want to introduce a dependency on a unit test framework into
    // this repo. Test 1: If x(0) = y(0) = 0, then x(t) = y(t) = 0 and z(t) =
    // z(0)exp(-βt).
    Real sigma = 10;
    Real beta = Real(8) / Real(3);
    Real rho = 28;
    Real tmax = 10;
    Real absolute_error = 1e-5;
    const std::array<Real, 3> initial_conditions{0, 0, Real(1)};
    auto solution = lorenz<double>(sigma, beta, rho, initial_conditions, tmax, absolute_error);
    auto const &skeleton = solution.states();

    for (auto const &s : skeleton)
    {
        Real t = s[0];
        Real x = s[1];
        Real y = s[2];
        Real z = s[3];
        if (abs(x) > std::numeric_limits<Real>::epsilon())
        {
            throw std::logic_error("x < eps doesn't hold");
        }
        if (abs(y) > std::numeric_limits<Real>::epsilon())
        {
            throw std::logic_error("y < eps doesn't hold");
        }
        Real expected = std::exp(-beta * t);
        if (abs(expected - z) > 100 * absolute_error)
        {
            std::cerr << "Expected z = " << expected << "\n";
            std::cerr << "Computed z = " << z << "\n";
            throw std::logic_error("z(t) = exp(-βt) doesn't hold");
        }
    }
}

#endif
