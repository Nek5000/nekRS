#include <iostream>
#include <vector>

#include <adios2.h>

int main(int argc, char *argv[])
{
    unsigned long N = 256;
    if (argc != 2)
    {
        std::cout << "Must pass adios config file!\n";
        return EXIT_FAILURE;
    }
    adios2::ADIOS adios(argv[1]);

    adios2::IO io = adios.DeclareIO("PluginTest");
    io.SetEngine("plugin");

    auto u = io.DefineVariable<double>("density", {N, N, N}, {0, 0, 0}, {N, N, N});
    adios2::Engine writer = io.Open("writer", adios2::Mode::Write);
    for (int64_t timeStep = 0; timeStep < 2; ++timeStep)
    {
        writer.BeginStep();
        std::vector<double> v(N * N * N, 3.2);
        std::cout << "Putting data at address   " << v.data() << " into inline writer.\n";
        writer.Put(u, v.data());
        writer.EndStep();
    }
    writer.Close();
    return EXIT_SUCCESS;
}
