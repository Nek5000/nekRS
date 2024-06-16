#if !defined(nekrs_adios_hpp_)
#define nekrs_adios_hpp_

#if !defined(NEKRS_ENABLE_ADIOS)

#error "No Adios2 installation was found"

#else

#define NEKRS_ADIOS_ENABLED
#include "platform.hpp"
#include "adios2.h"

class NekAdios
{
private:
  adios2::ADIOS *adios;
  adios2::IO io;
  adios2::Engine engine;

  using field = std::tuple<std::string, occa::memory, mesh_t *, dlong>;
  std::vector <field> userFieldList;

  std::vector<occa::memory> fldDataOut;
  occa::memory verticesOut;
  std::vector<occa::memory> work;

  const uint32_t VTK_CELL_TYPE = 12; // VTK_HEXAHEDRON

  bool uniform;

  bool movingMesh;

  bool initialized;
  std::string timerPrefix;

  mesh_t *mesh;
  mesh_t *mesh_vis;

  std::string vtkSchema()
  {
    std::string schema = R"(
    <VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
        <UnstructuredGrid>
            <Piece NumberOfPoints="NumOfPoints" NumberOfCells="NumOfCells">
                <Points>
                    <DataArray Name="vertices" />
                </Points>
                <Cells>
                    <DataArray Name="connectivity" />
                    <DataArray Name="types" />
                </Cells>
                <PointData>
    )";


    for (auto &entry : userFieldList) {
      const auto fieldName = std::get<0>(entry);
      schema += " <DataArray Name=\"" + fieldName + "\"/>\n";
    }

    schema += R"( <DataArray Name="TIME"> TIME </DataArray> )";

    schema += R"( 
                </PointData>
            </Piece>
        </UnstructuredGrid>
    </VTKFile>
  )";

    return schema;
  }

  template <typename T> adios2::Variable<T> defineVariable(const std::string &name, adios2::Dims dim = {})
  {
    auto var = io.InquireVariable<T>(name);
    if (var) {
      return var;
    }

    if (dim.size()) {
      return io.DefineVariable<T>(name, {}, {}, dim, adios2::ConstantDims);
    } else {
      return io.DefineVariable<T>(name);
    }
  }

public:
  NekAdios(){};

  NekAdios(mesh_t *meshIn,
           const std::string &name,
           const int N,
           const bool uniform_ = false,
           const std::string streamName = "default",
           const std::string configFile = "")
  {
    mesh = meshIn;
    uniform = uniform_;

    mesh_vis = [&]()
    {
      auto meshNew = meshIn;
      if (uniform || (N > 0 && N != mesh->N)) {
        meshNew = new mesh_t();
        meshNew->Nelements = mesh->Nelements;
        meshNew->dim = mesh->dim;
        meshNew->Nverts = mesh->Nverts;
        meshNew->Nfaces = mesh->Nfaces;
        meshNew->NfaceVertices = mesh->NfaceVertices;
        meshLoadReferenceNodesHex3D(meshNew, N, 0);
      }
      return meshNew;
    }();

    movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");

    timerPrefix = "nekAdios_" + streamName + "::";

    platform->timer.addUserStat(timerPrefix);

    if (configFile.size()) {
      adios = new adios2::ADIOS(configFile, platform->comm.mpiComm);
    } else {
      adios = new adios2::ADIOS(platform->comm.mpiComm);
    }
    io = adios->DeclareIO(streamName);
    io.SetEngine("BP5");
    io.DefineAttribute<uint32_t>("dimension", static_cast<uint32_t>(mesh_vis->dim));

    const auto fileName = name + ".bp";
    engine = io.Open(fileName, adios2::Mode::Write);

    initialized = true;
  }

  NekAdios(mesh_t *mesh,
           const std::string &name,
           const std::string streamName = "default",
           const std::string configFile = "") : 
           NekAdios(mesh, name, mesh->N, false, streamName, configFile) {}

  void close()
  {
    engine.Close();
  }

  ~NekAdios()
  {
    if (mesh_vis != mesh) {
      meshFree(mesh_vis);
    }
    close();
  }

  void addScalarField(const std::string& name, occa::memory o_fld, mesh_t *mesh_fld)
  {
    userFieldList.push_back(std::tuple{name, o_fld.slice(0, mesh_fld->Nlocal), mesh_fld, 0});
  }

  void addVectorField(const std::string& name, occa::memory o_fld, mesh_t *mesh_fld, dlong offset)
  {
    userFieldList.push_back(std::tuple{name, o_fld.slice(0, mesh_fld->dim * offset), mesh_fld, offset});
  }

  void clearFieldData()
  {
    userFieldList.clear();
  }

  template <typename OutputType>
  void convertField(mesh_t *mesh_fld, const std::vector<occa::memory> &o_fld, occa::memory& fldData)
  {
    const auto dim_fld = o_fld.size(); 
    nekrsCheck(dim_fld != 1 && dim_fld != 3,
               MPI_COMM_SELF,
               EXIT_FAILURE,
               "%s %d\n",
               "invalid field dimension=",
               dim_fld);

    auto fldDataPtr = static_cast<OutputType*>(fldData.ptr());
    for (int n = 0; n < fldData.size(); ++n) fldDataPtr[n] = 0;

    auto fld = [&](const int dim_i) 
    {
      nekrsCheck(o_fld[dim_i].dtype().name() != dfloatString,
                 MPI_COMM_SELF,
                 EXIT_FAILURE,
                 "%s %s\n",
                 "invalid field data type=",
                 o_fld[dim_i].dtype().name().c_str());

      if (!uniform && (mesh_vis->N == mesh_fld->N)) {
        return o_fld[dim_i]; 
      } 

      auto o_out = platform->o_memPool.reserve<dfloat>(mesh_vis->Nlocal);

      if (uniform) {
        mesh_fld->map2Uniform(mesh_vis, o_fld[dim_i], o_out);
      } else {
        mesh_fld->interpolate(mesh_vis, o_fld[dim_i], o_out);
      }

      o_out.copyTo(work[dim_i]);
      
      return work[dim_i];
    };

    if (dim_fld == 1) {
      auto fld0 = fld(0);
      auto fld0Ptr = static_cast<dfloat*>(fld0.ptr());

      for (int n = 0; n < fld0.size(); ++n) {
        fldDataPtr[n] = fld0Ptr[n];
      }
    } else if (dim_fld == 3) {
      auto fld0 = fld(0);
      auto fld0Ptr = static_cast<dfloat*>(fld0.ptr());
   
      auto fld1 = fld(1);
      auto fld1Ptr = static_cast<dfloat*>(fld1.ptr());

      auto fld2 = fld(2);
      auto fld2Ptr = static_cast<dfloat*>(fld2.ptr());

      // VTK expects AOS
      for (int n = 0; n < fld0.size(); ++n) {
        fldDataPtr[n * 3 + 0] = fld0Ptr[n];
        fldDataPtr[n * 3 + 1] = fld1Ptr[n];
        fldDataPtr[n * 3 + 2] = fld2Ptr[n];
      }
    }
  }

  template <typename OutputType = float>
  void write(double time, int tstep)
  {

    nekrsCheck(!initialized, MPI_COMM_SELF, EXIT_FAILURE, "%s\n", "called prior to initialization!");

    platform->timer.tic(timerPrefix + "write");

    if (platform->comm.mpiRank == 0) {
      std::cout << timerPrefix << " writing to " << engine.Name() << " ...";
    }
    const double tStart = MPI_Wtime();

    static bool firstTime = true;

    // pre-allocate all required host pool memory as 
    // no resizing is allowed prior to PerformPuts()
    verticesOut = platform->memPool.reserve<OutputType>(mesh_vis->dim * mesh_vis->Nlocal);
    for (auto &entry : userFieldList) {
      auto& o_fld = std::get<1>(entry);
      auto& offset_fld = std::get<3>(entry);
      const auto dim_fld = (offset_fld) ? o_fld.size() / offset_fld : 1;

      fldDataOut.push_back(platform->memPool.reserve<OutputType>(dim_fld * mesh_vis->Nlocal)); 
    }
    for (int i = 0; i < mesh_vis->dim; i++) {
      work.push_back(platform->memPool.reserve<dfloat>(mesh_vis->Nlocal));
    }

    if (firstTime) { 
      io.DefineAttribute<std::string>("vtk.xml", vtkSchema());
    }

    engine.BeginStep();

    size_t writtenBytes;

    auto var_time = defineVariable<double>("TIME");
    engine.Put(var_time, time);

    if (firstTime || movingMesh) {
      auto putMode = adios2::Mode::Sync;

      if (firstTime || movingMesh) {
        auto var_types = defineVariable<uint32_t>("types");
        engine.Put(var_types, VTK_CELL_TYPE, putMode);

        const uint32_t NumOfCells = mesh_vis->Nelements * std::pow(mesh_vis->N, mesh_vis->dim);
        auto var_NumOfCells = defineVariable<uint32_t>("NumOfCells");
        engine.Put(var_NumOfCells, NumOfCells, putMode);

        auto var_NumOfPoints = defineVariable<uint32_t>("NumOfPoints");
        engine.Put(var_NumOfPoints, static_cast<uint32_t>(mesh_vis->Nlocal), putMode);

        auto var_connectivity = defineVariable<uint64_t>("connectivity",
            {static_cast<size_t>(NumOfCells), static_cast<size_t>(mesh_vis->Nverts + 1)});

        const auto etov = [&] ()
        {
          std::vector<uint64_t> etov(NumOfCells * (mesh_vis->Nverts + 1));
          const auto Nq = mesh_vis->Nq;
          const auto NqPlane = Nq * Nq;

          auto it = etov.begin();
          for (int e = 0; e < mesh_vis->Nelements; ++e) {
            for (int z = 0; z < mesh_vis->N; ++z) {
              for (int y = 0; y < mesh_vis->N; ++y) {
                for (int x = 0; x < mesh_vis->N; ++x) {
                  it[0] = mesh_vis->Nverts;

                  // VTK_HEXAHEDRON ordering
                  it[1] = ((e * Nq + z) * Nq + y) * Nq + x;
                  it[2] = it[1] + 1;
                  it[3] = it[1] + 1 + Nq;
                  it[4] = it[1] + Nq;

                  it[5] = it[1] + NqPlane;
                  it[6] = it[2] + NqPlane;
                  it[7] = it[3] + NqPlane;
                  it[8] = it[4] + NqPlane;

                  it += (mesh_vis->Nverts + 1);
                }
              }
            }
          }

          return etov;
        }();

        engine.Put(var_connectivity, etov.data(), putMode);
        writtenBytes += etov.size() * sizeof(uint64_t);
      }

      std::vector<occa::memory> o_xyz;
      o_xyz.push_back(mesh->o_x);
      o_xyz.push_back(mesh->o_y);
      o_xyz.push_back(mesh->o_z);

      convertField<OutputType>(mesh, o_xyz, verticesOut);

      auto var_vertices =
          defineVariable<OutputType>("vertices", {static_cast<size_t>(mesh_vis->Nlocal), static_cast<size_t>(mesh_vis->dim)});

      engine.Put(var_vertices, static_cast<OutputType*>(verticesOut.ptr()));
      writtenBytes += verticesOut.size() * sizeof(OutputType);

      firstTime = false; 
    }

    int fldIdx = 0;
    for (auto &entry : userFieldList) {
      auto& fieldName = std::get<0>(entry);
      auto& o_fld = std::get<1>(entry);
      auto& mesh_fld = std::get<2>(entry);
      auto& offset_fld = std::get<3>(entry);

      std::vector<occa::memory> o_fldVec;
      const auto dim_fld = (offset_fld) ? o_fld.size() / offset_fld : 1;
      for (int i = 0; i < dim_fld; i++) {
        o_fldVec.push_back(o_fld.slice(i * offset_fld, mesh_fld->Nlocal));
      }

      auto& fldDataOutEntry = fldDataOut[fldIdx++];
      convertField<OutputType>(mesh_fld, o_fldVec, fldDataOutEntry);
 
      const auto count = [&]() {
        size_t dim_fld = fldDataOutEntry.size() / mesh_vis->Nlocal;
        if (dim_fld > 1) {
          return adios2::Dims{static_cast<size_t>(mesh_vis->Nlocal), dim_fld};
        } else {
          return adios2::Dims{static_cast<size_t>(mesh_vis->Nlocal)};
        }
      }();

      auto var = defineVariable<OutputType>(fieldName, count);
      engine.Put(var, static_cast<OutputType*>(fldDataOutEntry.ptr()));
      writtenBytes += fldDataOutEntry.size() * sizeof(OutputType);
    }

    engine.PerformPuts();
    engine.EndStep();

    verticesOut.free();
    fldDataOut.clear();
    work.clear();

    platform->timer.toc(timerPrefix + "write");

    MPI_Barrier(platform->comm.mpiComm);
    if (platform->comm.mpiRank == 0) {
      const auto timeWrite = MPI_Wtime() - tStart;
      printf(" done (%gs, %.2fGB/s)\n", timeWrite, writtenBytes/timeWrite/1e9);
    }
    fflush(stdout);

    initialized = true;
  }


};

#endif
#endif
