### ADIOS2 Lorenz System Example

The [Lorenz system](https://en.wikipedia.org/wiki/Lorenz_system) is a very simple ODE which exhibits large sensitivity
to initial conditions.

This example exhibits a number of capabilities of ADIOS: For example, how do we write data that is completely
independent of other data rapidly?
How do we write complex data structures that are not naturally layed out as matrices?

The ODE to be solved is

![Alt text](./lorenz_ode.svg)

We will use Taylor methods for ODEs, which are very easy to understand.

Namely, if we differentiate the right-hand side we get

![Alt text](./second_derivative.svg)

and we can propagate our solution via analytic continuation

![Alt text](./taylor_series.svg)

Normally, we would store the list

![Alt text](./ode_tuples.svg)

but we will be able to store less total data if we store the phase-space state

![Alt text](./ode_phase_space.svg)

and use cubic Hermite interpolation for rendering.
In addition, it doesn't make much sense to *solve* our equation more accurately than we render it, so we'll simply use
the second derivative to pick stepsizes.

A natural C++ data structure for this is a `std::vector<std::array<Real, 7>>`.
This will give us a nontrivial data structure to examine in ADIOS2.

Each MPI rank solves the Lorenz system with a different initial condition, and there is no coordination between ranks.
