/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Write a single HEXA8 cell (8 points in space for one box)
 * and provide a fides schema for visualizing it in ParaView 5.12 or later
 *
 * Created on: Nov 30, 2022
 *      Author: pnorbert
 */

#include <adios2.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>

int main(int argc, char *argv[])
{
    const size_t NPOINTS = 8;
    const size_t NCELLS = 1;
    const size_t NSTEPS = 10;
    adios2::ADIOS adios;

    const double X[NPOINTS] = {0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0};
    const double Y[NPOINTS] = {0.0, 0.0, 2.0, 2.0, 0.0, 0.0, 2.0, 2.0};
    const double Z[NPOINTS] = {0.0, 0.0, 0.0, 0.0, 3.0, 3.0, 3.0, 3.0};

    const int64_t cell[NCELLS * NPOINTS] = {0, 1, 2, 3, 4, 5, 6, 7};

    adios2::IO io = adios.DeclareIO("Output");
    io.SetEngine("BPFile");

    /*
     * Define global array: type, name, global dimensions
     * The local process' part (start, count) can be defined now or later
     * before Write().
     */
    adios2::Variable<double> varX =
        io.DefineVariable<double>("Mesh/Coords/X", {NPOINTS}, {0}, {NPOINTS});
    adios2::Variable<double> varY =
        io.DefineVariable<double>("Mesh/Coords/Y", {NPOINTS}, {0}, {NPOINTS});
    adios2::Variable<double> varZ =
        io.DefineVariable<double>("Mesh/Coords/Z", {NPOINTS}, {0}, {NPOINTS});
    adios2::Variable<int64_t> varCells =
        io.DefineVariable<int64_t>("Mesh/Cells", {NCELLS * NPOINTS}, {0}, {NCELLS * NPOINTS});

    adios2::Variable<size_t> varStep = io.DefineVariable<size_t>("step");
    adios2::Variable<double> varTime = io.DefineVariable<double>("time");

    adios2::Variable<double> varPointData =
        io.DefineVariable<double>("DataOnPoints", {NPOINTS}, {0}, {NPOINTS});
    adios2::Variable<double> varCellData =
        io.DefineVariable<double>("DataOnCells", {NCELLS}, {0}, {NCELLS});

    io.DefineAttribute<size_t>("npoints", NPOINTS);
    io.DefineAttribute<size_t>("ncells", NCELLS);

    adios2::Engine writer = io.Open("onecell.bp", adios2::Mode::Write);

    double pointdata[NPOINTS] = {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};
    double celldata[NCELLS] = {0.0};

    for (size_t step = 0; step < NSTEPS; step++)
    {
        writer.BeginStep();

        for (size_t p = 0; p < NPOINTS; p++)
        {
            pointdata[p] = std::min(1.0, pointdata[p] + 0.1);
        }
        for (size_t c = 0; c < NCELLS; c++)
        {
            celldata[c] += 0.1;
        }

        if (!step)
        {
            writer.Put(varX, X);
            writer.Put(varY, Y);
            writer.Put(varZ, Z);
            writer.Put(varCells, cell);
        }

        writer.Put<double>(varPointData, pointdata);
        writer.Put<double>(varCellData, celldata);
        writer.Put<size_t>(varStep, step);
        writer.Put<double>(varTime, static_cast<double>(step));

        writer.EndStep();
    }

    writer.Close();
    return 0;
}
