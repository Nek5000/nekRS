/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */

#include <iostream>
#include <vector>

#include <adios2.h>

int main(int argc, char *argv[])
{
    unsigned long N = 256;
    adios2::ADIOS adios;
    adios2::IO io = adios.DeclareIO("InlineReadWrite");
    io.SetEngine("Inline");
    auto u = io.DefineVariable<double>("u", {}, {}, {N}, adios2::ConstantDims);
    adios2::Engine writer = io.Open("writer", adios2::Mode::Write);
    adios2::Engine reader = io.Open("reader", adios2::Mode::Read);
    for (int64_t timeStep = 0; timeStep < 2; ++timeStep)
    {
        writer.BeginStep();
        std::vector<double> v(N, 3.2);
        std::cout << "Putting data at address   " << v.data() << " into inline writer.\n";
        writer.Put(u, v.data());
        writer.EndStep();

        reader.BeginStep();
        double *data = nullptr;
        reader.Get(u, &data);
        std::cout << "Getting data from address " << data << " via inline reader\n";
        reader.EndStep();
    }
    writer.Close();
    reader.Close();
    return 0;
}
