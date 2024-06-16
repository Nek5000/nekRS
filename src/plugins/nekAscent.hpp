#if !defined(nekrs_nekascent_hpp_)
#define nekrs_nekascent_hpp_

#if !defined(NEKRS_ENABLE_ASCENT)

#error "No Ascent installation was found"

#else

#define NEKRS_ASCENT_ENABLED
#include "platform.hpp"
#include "ascent.hpp"
#include "vtkh/vtkh.hpp"

namespace nekAscent
{

void setup(mesh_t *mesh_,
           const std::string& inputFile,
           int Nin_ = 0,
           bool uniform_ = false);
void run(const double time, const int tstep);
void finalize();
void addScalarField(const std::string& name, occa::memory o_fld, mesh_t *mesh_fld);
void addVectorField(const std::string& name, occa::memory o_fld, mesh_t *mesh_fld, dlong offset);
void clearData();

ascent::Ascent mAscent;
} // namespace nekAscent

namespace
{
conduit::Node mesh_data;

using field = std::tuple<std::string, occa::memory, mesh_t *, dlong>;
std::vector<field> userFieldList;
occa::memory o_connectivity;

mesh_t *mesh_in;
mesh_t *mesh_vis;

bool setupCalled = false;
bool updateMesh = true;

bool interpolate = false;
bool uniform = false;

std::string actionFile;
} // namespace

static void errHandler(const std::string &msg,
                       const std::string &file,
                       int line)
{
  nekrsAbort(MPI_COMM_SELF, EXIT_FAILURE, 
             "%s\n", msg.c_str());
}

static void initializeAscent()
{
  const double tStart = MPI_Wtime();

  MPI_Comm comm;
  MPI_Comm_dup(platform->comm.mpiComm, &comm);

  conduit::utils::set_warning_handler(errHandler);
  conduit::utils::set_error_handler(errHandler);

  conduit::Node ascent_opts;
  ascent_opts["mpi_comm"] = MPI_Comm_c2f(comm);
  //ascent_opts["exceptions"] = "forward";
  //ascent_opts["messages"] = "verbose";

  nekAscent::mAscent.open(ascent_opts);

  if (platform->comm.mpiRank == 0) {
    conduit::Node about;
    ascent::about(about);
    about.remove_child("license");
    about.remove_child("annotations");
    about.remove_child("git_sha1_abbrev");
    about.remove_child("git_tag");
    about.remove_child("compilers");
    about.remove_child("platform");
    about.remove_child("system");
    about.remove_child("web_client_root");
    about.remove_child("default_runtime");
    about["runtimes/ascent"].remove_child("status");
    std::cout << "---------------- Ascent.about() ----------------";
    std::cout << about.to_yaml() << std::endl;

    std::cout << vtkh::AboutVTKH() << std::endl;
  }

  platform->timer.set("nekAscent::setup::initializeAscent", MPI_Wtime() - tStart);

  fflush(stdout);
}

static void updateFieldData(occa::memory& o_fields)
{
  platform->timer.tic("nekAscent::run::update");

  auto fieldOffsetScan = [&]()
  {
    std::vector<dlong> offsetScan(userFieldList.size() + 1);
    offsetScan[0] = 0;

    int ifld = 1;
    for (auto &entry : userFieldList) {
      auto fieldName = std::get<0>(entry);
      auto o_fld = std::get<1>(entry);
      auto mesh_fld = std::get<2>(entry);
      auto& offset_fld = std::get<3>(entry);
      const auto dim_fld = (offset_fld) ? o_fld.size() / offset_fld : 1;
 
      offsetScan[ifld] =
          offsetScan[ifld - 1] +
          alignStride<dfloat>(mesh_vis->Np * (mesh_fld->Nelements + mesh_fld->totalHaloPairs)) * dim_fld;
 
      ifld++;
    }
    return offsetScan;
  }();

  if (interpolate) {
    if (!o_fields.isInitialized())
      o_fields = platform->o_memPool.reserve<dfloat>(fieldOffsetScan.back());
 
    const bool movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");
    if (updateMesh || movingMesh) {
      if (uniform) {
        mesh_in->map2Uniform(mesh_vis, mesh_in->o_x, mesh_vis->o_x);
        mesh_in->map2Uniform(mesh_vis, mesh_in->o_y, mesh_vis->o_y);
        mesh_in->map2Uniform(mesh_vis, mesh_in->o_z, mesh_vis->o_z);
      } else {
        mesh_in->interpolate(mesh_vis, mesh_in->o_x, mesh_vis->o_x);
        mesh_in->interpolate(mesh_vis, mesh_in->o_y, mesh_vis->o_y);
        mesh_in->interpolate(mesh_vis, mesh_in->o_z, mesh_vis->o_z);
      }
      updateMesh = false;
    }
  }

  int ifld = 0;
  const std::string str_uvw = "uvw";
  const std::string str_xyz = "xyz";
  for (auto &entry : userFieldList) {
    auto fieldName = std::get<0>(entry);
    auto o_fldIn = std::get<1>(entry);
    auto mesh_fld = std::get<2>(entry);
    auto& offset_fldIn = std::get<3>(entry);

    const auto dim_fld = (offset_fldIn) ? o_fldIn.size() / offset_fldIn : 1;
    const dlong offset_fld = (interpolate) ? fieldOffsetScan[ifld] / dim_fld : 0;
    const dlong fieldLength = mesh_fld->Nelements * mesh_vis->Np;

    auto o_fld = (interpolate) ? o_fields.slice(fieldOffsetScan[ifld]) : o_fldIn;

    if (dim_fld==1) { // scalar field
      if (interpolate) {
        if (uniform) {
          mesh_fld->map2Uniform(mesh_vis, o_fldIn, o_fld);
        } else {
          mesh_fld->interpolate(mesh_vis, o_fldIn, o_fld);
        }
      }
      mesh_data["fields/" + fieldName + "/association"] = "vertex";
      mesh_data["fields/" + fieldName + "/topology"] = "mesh";
      mesh_data["fields/" + fieldName + "/values"].set_external((dfloat *)o_fld.ptr(), fieldLength);

    } else if (dim_fld==2 || dim_fld==3) { // vector field

      for (int idim=0;idim<dim_fld;idim++) {
        auto o_fldInComp = o_fldIn.slice(idim*offset_fldIn);
        auto o_fldComp = (interpolate) ? o_fld.slice(idim*offset_fld) : o_fldInComp;

        if (interpolate) {
          if (uniform) {
            mesh_fld->map2Uniform(mesh_vis, o_fldInComp, o_fldComp);
          } else {
            mesh_fld->interpolate(mesh_vis, o_fldInComp, o_fldComp);
          }
        }

        // scalar field
        std::string fieldNameXYZ = fieldName + "_" + str_xyz.at(idim);
        mesh_data["fields/" + fieldNameXYZ + "/association"] = "vertex";
        mesh_data["fields/" + fieldNameXYZ + "/topology"] = "mesh";
        mesh_data["fields/" + fieldNameXYZ + "/values"].set_external((dfloat *)o_fldComp.ptr(), fieldLength);

        // vector component
        if (idim==0) {
          mesh_data["fields/" + fieldName + "/association"] = "vertex";
          mesh_data["fields/" + fieldName + "/topology"] = "mesh";
        }
        mesh_data["fields/" + fieldName + "/values/" + str_uvw.at(idim)].set_external((dfloat *)o_fldComp.ptr(), fieldLength);
      }
    } else {
      if (platform->comm.mpiRank == 0) {
        nekrsCheck(0, MPI_COMM_SELF, EXIT_FAILURE, "Unsupported vector dim: %d\n", dim_fld);
      }
    }

    ifld++;
  }

  platform->timer.toc("nekAscent::run::update");
}

void nekAscent::addScalarField(const std::string& name, occa::memory o_fld, mesh_t *mesh_fld)
{
  userFieldList.push_back(std::tuple{name, o_fld.slice(0, mesh_fld->Nlocal), mesh_fld, 0});
}

void nekAscent::addVectorField(const std::string& name, occa::memory o_fld, mesh_t *mesh_fld, dlong offset)
{
  userFieldList.push_back(std::tuple{name, o_fld.slice(0, mesh_fld->dim * offset), mesh_fld, offset});
}

void nekAscent::clearData()
{
  userFieldList.clear();
}

void nekAscent::setup(mesh_t *mesh_,
                      const std::string& inputFile,
                      int Nin_,
                      bool uniform_)
{
  const auto verbose = platform->options.compareArgs("VERBOSE", "TRUE") ? 1 : 0;
  mesh_in = mesh_;
  uniform = uniform_;
  const int Nin = (Nin_) ? Nin_ : mesh_in->N;

  platform->timer.addUserStat("nekAscent::");

  interpolate = (uniform || (Nin != mesh_in->N));

  const auto tStart = MPI_Wtime();
  if (platform->comm.mpiRank == 0) {
    printf("initializing nekAscent ");
    if (interpolate) {
      printf("(Nviz=%d", Nin);
      if (uniform)  printf(" +uniform"); 
      printf(") ...\n");
      fflush(stdout);
    }
  }

  if (platform->comm.mpiRank == 0) {
    nekrsCheck(!fs::exists(inputFile), MPI_COMM_SELF, EXIT_FAILURE, "Cannot find %s\n", inputFile.c_str());
  }
  actionFile = inputFile;

  initializeAscent();

  mesh_vis = [&]()
  {
    auto mesh = mesh_in;
    if (interpolate) {
      mesh = new mesh_t();
      mesh->Nelements = mesh_in->Nelements;
      mesh->dim = mesh_in->dim;
      mesh->Nverts = mesh_in->Nverts;
      mesh->Nfaces = mesh_in->Nfaces;
      mesh->NfaceVertices = mesh_in->NfaceVertices;
      meshLoadReferenceNodesHex3D(mesh, Nin, 0);
 
      mesh->o_x = platform->device.malloc<dfloat>(mesh->Nlocal);
      mesh->o_y = platform->device.malloc<dfloat>(mesh->Nlocal);
      mesh->o_z = platform->device.malloc<dfloat>(mesh->Nlocal);
    }
    return mesh;
  }();

  mesh_data["coordsets/coords/type"] = "explicit";
  mesh_data["coordsets/coords/values/x"].set_external((dfloat *)mesh_vis->o_x.ptr(), mesh_vis->Nlocal);
  mesh_data["coordsets/coords/values/y"].set_external((dfloat *)mesh_vis->o_y.ptr(), mesh_vis->Nlocal);
  mesh_data["coordsets/coords/values/z"].set_external((dfloat *)mesh_vis->o_z.ptr(), mesh_vis->Nlocal);

  mesh_data["topologies/mesh/type"] = "unstructured";
  mesh_data["topologies/mesh/coordset"] = "coords";
  mesh_data["topologies/mesh/elements/shape"] = "hex";

  o_connectivity = [&]()
  {
    const dlong Nverts = mesh_vis->Nelements * std::pow(mesh_vis->N, mesh_vis->dim) * mesh_vis->Nverts;
    std::vector<dlong> etov(Nverts);
    auto o_etov = platform->device.malloc<dlong>(etov.size());

    auto it = etov.begin();
    for (int e = 0; e < mesh_vis->Nelements; ++e) {
      for (int z = 0; z < mesh_vis->N; ++z) {
        for (int y = 0; y < mesh_vis->N; ++y) {
          for (int x = 0; x < mesh_vis->N; ++x) {
            const dlong Nq = mesh_vis->Nq;
            it[0] = ((e * Nq + z) * Nq + y) * Nq + x;
            it[1] = it[0] + 1;
            it[2] = it[0] + Nq + 1;
            it[3] = it[0] + Nq;
            it[4] = it[0] + Nq * Nq;
            it[5] = it[1] + Nq * Nq;
            it[6] = it[2] + Nq * Nq;
            it[7] = it[3] + Nq * Nq;
            it += mesh_vis->Nverts;
          }
        }
      }
    }
    o_etov.copyFrom(etov.data());
    return o_etov;
  }();
  mesh_data["topologies/mesh/elements/connectivity"].set_external((dlong *)o_connectivity.ptr(), o_connectivity.size());

  const auto tSetup = MPI_Wtime() - tStart;
  platform->timer.set("nekAscent::setup", tSetup);
  if (platform->comm.mpiRank == 0) {
    printf("\ndone (%gs)\n\n", tSetup);
    if (verbose) {
      mesh_data.print();
    }
  }
  fflush(stdout);

  setupCalled = true;
}


void nekAscent::run(const double time, const int tstep)
{
  nekrsCheck(!setupCalled, MPI_COMM_SELF, EXIT_FAILURE, "%s\n", "called prior to nekAscent::setup()!");
  const auto verbose = platform->options.compareArgs("VERBOSE", "TRUE") ? 1 : 0;

  platform->timer.tic("nekAscent::run");
  if (platform->comm.mpiRank == 0) {
    std::cout << "processing " << actionFile << std::endl;
  }

  mesh_data["state/cycle"] = tstep;
  mesh_data["state/time"] = time;

  occa::memory o_work;
  updateFieldData(o_work);

  if (platform->comm.mpiRank == 0 && verbose) {
    std::cout << "---------------- Ascent mesh_data ----------------" << std::endl;
    mesh_data.print();
    fflush(stdout);
  }
  mAscent.publish(mesh_data);

  conduit::Node triggers;
  triggers["t1/params/condition"] = "True"; // control the condition in udf, not ascent
  triggers["t1/params/actions_file"] = actionFile;

  conduit::Node actions;
  conduit::Node &add_triggers = actions.append();
  add_triggers["action"] = "add_triggers";
  add_triggers["triggers"] = triggers;

  mAscent.execute(actions);

  o_work.free();

  platform->timer.toc("nekAscent::run");
}

void nekAscent::finalize()
{
  if (setupCalled) {
    o_connectivity.free();
    mAscent.close();
  }
}

#endif
#endif // hpp
