#include "nrs.hpp"
#include "nekrsAscent.hpp"
#include "platform.hpp"
#include "linAlg.hpp"

// private members
namespace {
nekrsAscent::fields userFieldList;
static mesh_t *mesh;
static int Nfields;
static dlong fieldOffset;
static dlong Nvertices;
static occa::memory o_fields;
static occa::memory o_Xcoord;
static occa::memory o_Ycoord;
static occa::memory o_Zcoord;
static occa::memory o_connectivity;
static std::vector<unsigned int> ghosts;
static bool setupCalled = false;
} // namespace


void nekrsAscent::setup(mesh_t *mesh_, const dlong fieldOffset_, const fields& flds) {
// TODO add flag for initialized

  const double tStart = MPI_Wtime();
  if (platform->comm.mpiRank == 0) {
    printf("Initialize Ascent ...\n");
    fflush(stdout);
  }

  mesh = mesh_;
  userFieldList = flds;
  Nfields = userFieldList.size();
  fieldOffset = fieldOffset_;

  const int verbose = platform->options.compareArgs("VERBOSE", "TRUE") ? 1 : 0;
  
  MPI_Comm comm;
  MPI_Comm_dup(platform->comm.mpiComm, &comm);

  conduit::Node ascent_opts;
  ascent_opts["mpi_comm"]=MPI_Comm_c2f(comm);

//  ascent_opts["runtime/type"] = "ascent";
/*
  std::string backend; // TODO: do nothing? openmp? dpcpp? support mismatch?
  platform->options.getArgs("THREAD MODEL", backend);
  if (backend == "CUDA") {
    ascent_opts["runtime/backend"] = "cuda";
  } else if (backend == "HIP") {
    ascent_opts["runtime/backend"] = "hip";
  } else if (backend == "CPU" || backend == "SERIAL") {
    ascent_opts["runtime/backend"] = "serial";
  } else {
    ascent_opts["runtime/backend"] = "openmp";
  }
*/
  if (verbose) { // FIXME: This doesn't do much?
    ascent_opts["ascent_info"] = "verbose";
    ascent_opts["messages"] = "verbose";
  }
  mAscent.open(ascent_opts);

  const double tSetup = MPI_Wtime() - tStart; 
  platform->timer.set("insituAscentSetup", tSetup);
  if (platform->comm.mpiRank == 0) {
    printf("done (%gs)\n\n", tSetup);
    std::cout << ascent::about() << std::endl; // TODO verbose?
  }
  fflush(stdout);

  dlong Ncells = mesh->Nelements * (mesh->Nq - 1) * (mesh->Nq - 1) * (mesh->Nq - 1);
  Nvertices = Ncells * 8;

  o_Xcoord = platform->device.malloc<dfloat>(fieldOffset);
  o_Ycoord = platform->device.malloc<dfloat>(fieldOffset);
  o_Zcoord = platform->device.malloc<dfloat>(fieldOffset);
  o_connectivity = platform->device.malloc<dlong>(Nvertices);
  o_fields = platform->device.malloc<dfloat>(Nfields*fieldOffset);

  // Copy Mesh
  o_Xcoord.copyFrom(mesh->o_x, mesh->Nlocal);
  o_Ycoord.copyFrom(mesh->o_y, mesh->Nlocal);
  o_Zcoord.copyFrom(mesh->o_z, mesh->Nlocal);

  // Calculate connectivity
  {
    std::vector<dlong> a_etov(Nvertices);
    auto it = a_etov.begin();
    for(int elem=0; elem<mesh->Nelements; ++elem)
      for(int z=0; z < mesh->Nq-1; ++z)
        for(int y=0; y < mesh->Nq-1; ++y)
          for(int x=0; x < mesh->Nq-1; ++x) {
                 it[0] = ((elem * mesh->Nq + z) * mesh->Nq + y) * mesh->Nq + x;
                 it[1] = it[0] + 1;
                 it[2] = it[0] + mesh->Nq+1;
                 it[3] = it[0] + mesh->Nq;
                 it[4] = it[0] + mesh->Nq*mesh->Nq;
                 it[5] = it[1] + mesh->Nq*mesh->Nq;
                 it[6] = it[2] + mesh->Nq*mesh->Nq;
                 it[7] = it[3] + mesh->Nq*mesh->Nq;
            it += 8;  
          }
    o_connectivity.copyFrom(a_etov.data(), a_etov.size());
  }

  // halo node FIXME: ??
  if( mesh->totalHaloPairs ) {
    ghosts.resize(mesh->Nlocal, 0);
    for (int n=0; n< mesh->totalHaloPairs * mesh->Nfp; n++) {
      ghosts[mesh->haloGetNodeIds[n]] = 1;
    }
  }

  setupCalled = true;
}

void nekrsAscent::run(const double time, const int tstep) { 

  nrsCheck(!setupCalled, MPI_COMM_SELF, EXIT_FAILURE,
           "%s\n", "called prior to nekrsAscent::setup()!");

  platform->timer.tic("insituAscentRun",1);

  conduit::Node mesh_data;
  mesh_data["state/cycle"] = tstep;
  mesh_data["state/time"] = time;

  // TODO: only copy mesh if moving mesh
  o_Xcoord.copyFrom(mesh->o_x, mesh->Nlocal);
  o_Ycoord.copyFrom(mesh->o_y, mesh->Nlocal);
  o_Zcoord.copyFrom(mesh->o_z, mesh->Nlocal);

  //coordinate system data
  mesh_data["coordsets/coords/type"] = "explicit";
  mesh_data["coordsets/coords/values/x"].set_external((dfloat*)o_Xcoord.ptr(), mesh->Nlocal);
  mesh_data["coordsets/coords/values/y"].set_external((dfloat*)o_Ycoord.ptr(), mesh->Nlocal);
  mesh_data["coordsets/coords/values/z"].set_external((dfloat*)o_Zcoord.ptr(), mesh->Nlocal);

  // topology data
  mesh_data["topologies/mesh/type"]           = "unstructured";
  mesh_data["topologies/mesh/coordset"]       = "coords";
  mesh_data["topologies/mesh/elements/shape"] = "hex";  // Note "hexs" - documentation on Ascent/Conduit webpage is incorrect TODO: double check this
  mesh_data["topologies/mesh/elements/connectivity"].set_external((dlong*) o_connectivity.ptr(), Nvertices);

  // halo nodes FIXME??
  if( mesh->totalHaloPairs ) {
    mesh_data["fields/ghosts/association"] = "vertex";
    mesh_data["fields/ghosts/topology"] = "mesh";
    mesh_data["fields/ghosts/values"].set(ghosts);
  }
  
  // fields
  platform->linAlg->fill(fieldOffset * Nfields, 0.0, o_fields);
  int ifld = 0;
  for(auto& entry : userFieldList) {
    std::string fieldName = std::get<0>(entry);
    occa::memory o_fieldValue = std::get<1>(entry);
    dlong fieldLength = std::get<2>(entry);

    auto o_fld = o_fields.slice(ifld*fieldOffset);
    o_fld.copyFrom(o_fieldValue, fieldLength);

    mesh_data["fields/" + fieldName + "/association"]  = "vertex";
    mesh_data["fields/" + fieldName + "/topology"]     = "mesh";
    mesh_data["fields/" + fieldName + "/values"].set_external((dfloat*) o_fieldValue.ptr(), fieldLength);
    ifld++;
  }

  mAscent.publish(mesh_data);
  conduit::Node actions;
  mAscent.execute(actions);

  platform->timer.toc("insituAscentRun"); 
  // TODO, how to print this? nekrsAscent::stat?
}

void nekrsAscent::printStat(const int tstep) { //TODO: try to extract img info??

  const int verbose = platform->options.compareArgs("VERBOSE", "TRUE") ? 1 : 0;
  int freq = 500, numSteps = -1;
  platform->options.getArgs("RUNTIME STATISTICS FREQUENCY", freq);
  platform->options.getArgs("NUMBER TIMESTEPS", numSteps);

  std::string tag = "insituAscentRun";
  std::string name = "    insituAscentRun     ";
  std::string type = "DEVICE:MAX";

  if (freq && tstep>0) {
    if (tstep % freq ==0 || tstep==numSteps) {
      const long long int nCalls = platform->timer.count(tag);
      const double tTag = platform->timer.query(tag, type);
      if (tTag > 0) {
        if (platform->comm.mpiRank == 0) {
          std::cout << name << tTag << "s";
          std::cout << "  " << nCalls << " calls\n";
        }
      }
    }
  }

  if (verbose) {
    if (platform->comm.mpiRank==0) { 
      conduit::Node ascent_info;
      mAscent.info(ascent_info);
      ascent_info.print();
          
      conduit::Node opts;
      opts["num_children_threshold"] = -1;
      opts["num_elements_threshold"] = -1;
      opts["depth"] = 1; 
      ascent_info.to_summary_string_stream(std::cout,opts);
      std::cout << std::endl;
    }
  }            
}


void nekrsAscent::finalize() {
  mAscent.close();
}
