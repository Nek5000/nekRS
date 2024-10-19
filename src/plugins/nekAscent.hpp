#if !defined(NEKRS_ENABLE_ASCENT)
#error "No Ascent installation was found"
#endif

#if !defined(nekrs_nekascent_hpp_)
#define nekrs_nekascent_hpp_

#define NEKRS_ASCENT_ENABLED
#include "platform.hpp"
#include "ascent.hpp"
#include "vtkh/vtkh.hpp"
#include "mesh.h"
#include <future>

namespace nekAscent
{

namespace
{
conduit::Node ascent_opts;
ascent::Ascent mAscent;
conduit::Node mesh_data;
conduit::Node triggers;
conduit::Node actions;

occa::memory o_connectivity;
occa::memory o_x, o_y, o_z;
occa::memory o_work;

mesh_t *mesh_in;
mesh_t *mesh_vis;

using field = std::tuple<std::string, std::vector<occa::memory>, mesh_t *>;
std::vector<field> userFieldList;

std::future<void> asyncRunner;

bool async = false;
bool stageThroughHost = false;
bool setupCalled = false;
bool interpolate = false;
bool uniform = false;

void errHandler(const std::string &msg, const std::string &file, int line)
{
  nekrsAbort(MPI_COMM_SELF, EXIT_FAILURE, "%s\n", msg.c_str());
}

void initializeAscent()
{
  const double tStart = MPI_Wtime();

  MPI_Comm comm;
  MPI_Comm_dup(platform->comm.mpiComm, &comm);

  conduit::utils::set_warning_handler(errHandler);
  conduit::utils::set_error_handler(errHandler);

  ascent_opts["mpi_comm"] = MPI_Comm_c2f(comm);
  // ascent_opts["runtime/vtkm/backend"] = "serial";
  //  ascent_opts["exceptions"] = "forward";
  //  ascent_opts["messages"] = "verbose";

  mAscent.open(ascent_opts);

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

  fflush(stdout);
}

void allocWork(size_t size)
{
  if (!o_work.isInitialized()) {
    // cannot use memPool as potential realloc can lead to dangling points when executing in async mode
    if (stageThroughHost) {
      o_work = platform->device.mallocHost<dfloat>(size);
    } else {
      o_work = platform->device.malloc<dfloat>(size);
    }
  } else {
    nekrsCheck(o_work.size() < size,
               MPI_COMM_SELF,
               EXIT_FAILURE,
               "%s\n",
               "Invalid attempt to realloc o_work!");
  }
}

void updateFieldData()
{
  static auto firstTime = true;

  platform->timer.tic("nekAscent::run::update");

  const auto fieldOffsetScan = [&]() {
    std::vector<size_t> offsetScan(userFieldList.size() + 1);

    int ifld = 1;
    offsetScan[0] = 0;
    for (auto &entry : userFieldList) {
      auto o_fld = std::get<1>(entry);
      const auto dim_fld = o_fld.size();

      offsetScan[ifld] =
          offsetScan[ifld - 1] + alignStride<dfloat>(mesh_vis->Np * mesh_vis->Nelements) * dim_fld;

      ifld++;
    }
    return offsetScan;
  }();

  const auto movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");
  if (firstTime || movingMesh) {
    if (interpolate || uniform) {
      if (uniform) {
        mesh_in->interpolate(mesh_in->o_x, mesh_vis, mesh_vis->o_x, true);
        mesh_in->interpolate(mesh_in->o_y, mesh_vis, mesh_vis->o_y, true);
        mesh_in->interpolate(mesh_in->o_z, mesh_vis, mesh_vis->o_z, true);
      } else {
        mesh_in->interpolate(mesh_in->o_x, mesh_vis, mesh_vis->o_x);
        mesh_in->interpolate(mesh_in->o_y, mesh_vis, mesh_vis->o_y);
        mesh_in->interpolate(mesh_in->o_z, mesh_vis, mesh_vis->o_z);
      }
    }

    if (stageThroughHost || async) {
      o_x.copyFrom(mesh_vis->o_x);
      o_y.copyFrom(mesh_vis->o_y);
      o_z.copyFrom(mesh_vis->o_z);
    }
  }

  int ifld = 0;
  for (auto &entry : userFieldList) {
    const auto &fieldName = std::get<0>(entry);
    const auto &o_fldIn = std::get<1>(entry);
    const auto &mesh_fld = std::get<2>(entry);

    const auto dim_fld = o_fldIn.size();
    const auto fieldOffsetWork = (fieldOffsetScan[ifld + 1] - fieldOffsetScan[ifld]) / dim_fld;

    for (int idim = 0; idim < dim_fld; idim++) {
      auto data = [&]() {
        occa::memory o_fldOut;
        if (interpolate || uniform) {
          auto o_tmp = platform->deviceMemoryPool.reserve<dfloat>(mesh_vis->Nlocal);
          if (uniform) {
            mesh_fld->interpolate(o_fldIn.at(idim), mesh_vis, o_tmp, true);
          } else {
            mesh_fld->interpolate(o_fldIn.at(idim), mesh_vis, o_tmp);
          }
          allocWork(fieldOffsetScan.back());
          auto o_fldWork = o_work.slice(fieldOffsetScan[ifld]);
          o_fldOut = o_fldWork.slice(idim * fieldOffsetWork, o_tmp.size());
          o_fldOut.copyFrom(o_tmp);
        } else if (async) {
          allocWork(fieldOffsetScan.back());
          auto o_fldWork = o_work.slice(fieldOffsetScan[ifld]);
          o_fldOut = o_fldWork.slice(idim * fieldOffsetWork, mesh_vis->Nlocal);
          o_fldOut.copyFrom(o_fldIn.at(idim), o_fldIn.at(idim).size());
        } else {
          occa::memory o_tmp;
          if (stageThroughHost) {
            o_tmp = platform->memoryPool.reserve<dfloat>(mesh_vis->Nlocal);
            auto ptr = o_tmp.ptr<dfloat>();
            for (int i = 0; i < o_tmp.size(); i++) {
              ptr[i] = 0.0;
            }
          } else {
            o_tmp = platform->deviceMemoryPool.reserve<dfloat>(mesh_vis->Nlocal);
            platform->linAlg->fill(o_tmp.size(), 0.0, o_tmp);
          }
          o_tmp.copyFrom(o_fldIn.at(idim), o_fldIn.at(idim).size());
          o_fldOut = o_tmp;
        }

        return o_fldOut;
      }();

      {
        const std::string str_xyz = "xyz";
// workaround for https://github.com/Alpine-DAV/ascent/issues/1329
#if 1
        const auto fieldNameXYZ = (dim_fld > 1) ? fieldName + "_" + str_xyz.at(idim) : fieldName;
        mesh_data["fields/" + fieldNameXYZ + "/association"] = "vertex";
        mesh_data["fields/" + fieldNameXYZ + "/topology"] = "mesh";
        mesh_data["fields/" + fieldNameXYZ + "/values"].set_external(data.ptr<dfloat>(), data.size());
#else
        const auto fieldNameXYZ = fieldName;
        if (idim == 0) {
          mesh_data["fields/" + fieldNameXYZ + "/association"] = "vertex";
          mesh_data["fields/" + fieldNameXYZ + "/topology"] = "mesh";
        }
        mesh_data["fields/" + fieldNameXYZ + "/values/" + str_xyz.at(idim)].set_external(data.ptr<dfloat>(),
                                                                                         data.size());
#endif
      }
    }

    ifld++;
  }

  mAscent.publish(mesh_data);

  platform->timer.toc("nekAscent::run::update");

  if (platform->comm.mpiRank == 0 && platform->verbose) {
    std::cout << "---------------- Ascent mesh_data ----------------" << std::endl;
    conduit::Node mesh_copy;
    mesh_copy.set(mesh_data);
    mesh_copy.print();
    fflush(stdout);
  }

  firstTime = false;
}

} // namespace

void addVariable(const std::string &name, mesh_t *mesh_fld, const std::vector<deviceMemory<dfloat>> &fld)
{
  std::vector<occa::memory> fld_;
  for (const auto &entry : fld) {
    fld_.push_back(entry);
  }

  const auto dim = fld_.size();
  nekrsCheck(dim != 1 && dim != 3, MPI_COMM_SELF, EXIT_FAILURE, "Unsupported vector dim: %d\n", dim);

  userFieldList.push_back(std::tuple{name, fld_, mesh_fld});
}

void clearData()
{
  userFieldList.clear();
}

void setup(mesh_t *mesh_,
           const std::string &actionFile,
           int Nin_ = 0,
           bool uniform_ = false,
           bool stageThroughHost_ = false bool async_ = false)
{
  mesh_in = mesh_;
  const int Nin = (Nin_) ? Nin_ : mesh_in->N;
  uniform = uniform_;
  async = async_;
  stageThroughHost = stageThroughHost_;
  if (platform->serial) {
    stageThroughHost = false;
  }

  if (stageThroughHost) {
    async = true;
  }

  if (async) {
    int provided;
    MPI_Query_thread(&provided);
    nekrsCheck(provided != MPI_THREAD_MULTIPLE,
               MPI_COMM_SELF,
               EXIT_FAILURE,
               "%s\n",
               "running Ascent in async mode requires NEKRS_MPI_THREAD_MULTIPLE=1");
  }

  platform->timer.addUserStat("nekAscent::");

  interpolate = (Nin != mesh_in->N);

  const auto tStart = MPI_Wtime();
  if (platform->comm.mpiRank == 0) {
    printf("initializing nekAscent ");
    if (interpolate || uniform) {
      printf("(Nviz=%d", Nin);
      if (uniform) {
        printf(" +uniform");
      }
      printf(") ...\n");
      fflush(stdout);
    }
  }

  if (platform->comm.mpiRank == 0) {
    nekrsCheck(!fs::exists(actionFile), MPI_COMM_SELF, EXIT_FAILURE, "Cannot find %s\n", actionFile.c_str());
  }

  initializeAscent();

  mesh_vis = [&]() {
    auto mesh = mesh_in;
    o_x = mesh->o_x;
    o_y = mesh->o_y;
    o_z = mesh->o_z;

    if (interpolate || uniform) {
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
      o_x = mesh->o_x;
      o_y = mesh->o_y;
      o_z = mesh->o_z;
    } else if (async && !stageThroughHost) {
      o_x = platform->device.malloc<dfloat>(mesh->Nlocal);
      o_y = platform->device.malloc<dfloat>(mesh->Nlocal);
      o_z = platform->device.malloc<dfloat>(mesh->Nlocal);
    }

    if (stageThroughHost) {
      o_x = platform->device.mallocHost<dfloat>(mesh->Nlocal);
      o_y = platform->device.mallocHost<dfloat>(mesh->Nlocal);
      o_z = platform->device.mallocHost<dfloat>(mesh->Nlocal);
    }

    return mesh;
  }();

  o_connectivity = [&]() {
    const dlong Nverts = mesh_vis->Nelements * std::pow(mesh_vis->N, mesh_vis->dim) * mesh_vis->Nverts;
    std::vector<dlong> etov(Nverts);

    occa::memory o_etov;
    if (stageThroughHost) {
      o_etov = platform->device.mallocHost<dlong>(etov.size());
    } else {
      o_etov = platform->device.malloc<dlong>(etov.size());
    }

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

  mesh_data["coordsets/coords/type"] = "explicit";
  mesh_data["coordsets/coords/values/x"].set_external(o_x.ptr<dfloat>(), o_x.size());
  mesh_data["coordsets/coords/values/y"].set_external(o_y.ptr<dfloat>(), o_y.size());
  mesh_data["coordsets/coords/values/z"].set_external(o_z.ptr<dfloat>(), o_z.size());

  mesh_data["topologies/mesh/type"] = "unstructured";
  mesh_data["topologies/mesh/coordset"] = "coords";
  mesh_data["topologies/mesh/elements/shape"] = "hex";
  mesh_data["topologies/mesh/elements/connectivity"].set_external((dlong *)o_connectivity.ptr(),
                                                                  o_connectivity.size());

  triggers["t1/params/condition"] = "True";
  triggers["t1/params/actions_file"] = actionFile;

  conduit::Node &add_triggers = actions.append();
  add_triggers["action"] = "add_triggers";
  add_triggers["triggers"] = triggers;

  const auto tSetup = MPI_Wtime() - tStart;
  platform->timer.set("nekAscent::setup", tSetup);
  if (platform->comm.mpiRank == 0) {
    printf("\ndone (%gs)\n\n", tSetup);
    if (platform->verbose) {
      conduit::Node mesh_copy;
      mesh_copy.set(mesh_data);
      mesh_copy.print();
    }
  }
  fflush(stdout);

  setupCalled = true;
}

void run(const double time, const int tstep)
{
  nekrsCheck(!setupCalled, MPI_COMM_SELF, EXIT_FAILURE, "%s\n", "called prior to nekAscent::setup()!");

  platform->timer.tic("nekAscent::run");

  if (asyncRunner.valid()) {
    platform->timer.tic("nekAscent::run::execute");
    asyncRunner.wait(); // wait until previous task is completed
    platform->timer.toc("nekAscent::run::execute");
  }

  if (platform->comm.mpiRank == 0) {
    std::cout << "processing Ascent action file ..." << std::endl << std::flush;
  }

  mesh_data["state/cycle"] = tstep;
  mesh_data["state/time"] = time;

  updateFieldData();

  if (async) {
    asyncRunner = std::async(std::launch::async, [&]() { mAscent.execute(actions); });
  } else {
    platform->timer.tic("nekAscent::run::execute");
    mAscent.execute(actions);
    platform->timer.toc("nekAscent::run::execute");
  }

  platform->timer.toc("nekAscent::run");
}

void finalize()
{
  if (asyncRunner.valid()) {
    platform->timer.tic("nekAscent::run::execute");
    asyncRunner.get();
    platform->timer.toc("nekAscent::run::execute");
  }

  if (setupCalled) {
    o_connectivity.free();
    mAscent.close();
  }

  o_work.free();
}

} // namespace nekAscent

#endif // hpp
