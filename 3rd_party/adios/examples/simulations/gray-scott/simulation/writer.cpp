/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */

#include "writer.h"

void define_bpvtk_attribute(const Settings &s, adios2::IO &io)
{
    auto lf_VTKImage = [](const Settings &s, adios2::IO &io) {
        const std::string extent = "0 " + std::to_string(s.L) + " " + "0 " + std::to_string(s.L) +
                                   " " + "0 " + std::to_string(s.L);

        const std::string imageData = R"(
        <?xml version="1.0"?>
        <VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">
          <ImageData WholeExtent=")" + extent +
                                      R"(" Origin="0 0 0" Spacing="1 1 1">
            <Piece Extent=")" + extent +
                                      R"(">
              <CellData Scalars="U">
                  <DataArray Name="U" />
                  <DataArray Name="V" />
                  <DataArray Name="TIME">
                    step
                  </DataArray>
              </CellData>
            </Piece>
          </ImageData>
        </VTKFile>)";

        io.DefineAttribute<std::string>("vtk.xml", imageData);
    };

    if (s.mesh_type == "image")
    {
        lf_VTKImage(s, io);
    }
    else if (s.mesh_type == "structured")
    {
        throw std::invalid_argument("ERROR: mesh_type=structured not yet "
                                    "   supported in settings.json, use mesh_type=image instead\n");
    }
    // TODO extend to other formats e.g. structured
}

Writer::Writer(const Settings &settings, const GrayScott &sim, adios2::IO io)
: settings(settings), io(io)
{
    io.DefineAttribute<double>("F", settings.F);
    io.DefineAttribute<double>("k", settings.k);
    io.DefineAttribute<double>("dt", settings.dt);
    io.DefineAttribute<double>("Du", settings.Du);
    io.DefineAttribute<double>("Dv", settings.Dv);
    io.DefineAttribute<double>("noise", settings.noise);
    // define VTK visualization schema as an attribute
    if (!settings.mesh_type.empty())
    {
        define_bpvtk_attribute(settings, io);
    }

    // add attributes for Fides
    io.DefineAttribute<std::string>("Fides_Data_Model", "uniform");
    double origin[3] = {0.0, 0.0, 0.0};
    io.DefineAttribute<double>("Fides_Origin", &origin[0], 3);
    double spacing[3] = {0.1, 0.1, 0.1};
    io.DefineAttribute<double>("Fides_Spacing", &spacing[0], 3);
    io.DefineAttribute<std::string>("Fides_Dimension_Variable", "U");

    std::vector<std::string> varList = {"U", "V"};
    std::vector<std::string> assocList = {"points", "points"};
    io.DefineAttribute<std::string>("Fides_Variable_List", varList.data(), varList.size());
    io.DefineAttribute<std::string>("Fides_Variable_Associations", assocList.data(),
                                    assocList.size());

    var_u = io.DefineVariable<double>("U", {settings.L, settings.L, settings.L},
                                      {sim.offset_z, sim.offset_y, sim.offset_x},
                                      {sim.size_z, sim.size_y, sim.size_x});

    var_v = io.DefineVariable<double>("V", {settings.L, settings.L, settings.L},
                                      {sim.offset_z, sim.offset_y, sim.offset_x},
                                      {sim.size_z, sim.size_y, sim.size_x});

    if (settings.adios_memory_selection)
    {
        var_u.SetMemorySelection({{1, 1, 1}, {sim.size_z + 2, sim.size_y + 2, sim.size_x + 2}});
        var_v.SetMemorySelection({{1, 1, 1}, {sim.size_z + 2, sim.size_y + 2, sim.size_x + 2}});
    }

    var_step = io.DefineVariable<int>("step");
}

void Writer::open(const std::string &fname, bool append)
{
    adios2::Mode mode = adios2::Mode::Write;
    if (append)
    {
        mode = adios2::Mode::Append;
    }
    writer = io.Open(fname, mode);
}

void Writer::write(int step, const GrayScott &sim)
{
    if (!sim.size_x || !sim.size_y || !sim.size_z)
    {
        writer.BeginStep();
        writer.EndStep();
        return;
    }

    if (settings.adios_memory_selection)
    {
        const std::vector<double> &u = sim.u_ghost();
        const std::vector<double> &v = sim.v_ghost();

        writer.BeginStep();
        writer.Put<int>(var_step, &step);
        writer.Put<double>(var_u, u.data());
        writer.Put<double>(var_v, v.data());
        writer.EndStep();
    }
    else if (settings.adios_span)
    {
        writer.BeginStep();

        writer.Put<int>(var_step, &step);

        // provide memory directly from adios buffer
        adios2::Variable<double>::Span u_span = writer.Put<double>(var_u);
        adios2::Variable<double>::Span v_span = writer.Put<double>(var_v);

        // populate spans
        sim.u_noghost(u_span.data());
        sim.v_noghost(v_span.data());

        writer.EndStep();
    }
    else
    {
        std::vector<double> u = sim.u_noghost();
        std::vector<double> v = sim.v_noghost();

        writer.BeginStep();
        writer.Put<int>(var_step, &step);
        writer.Put<double>(var_u, u.data());
        writer.Put<double>(var_v, v.data());
        writer.EndStep();
    }
}

void Writer::close() { writer.Close(); }
