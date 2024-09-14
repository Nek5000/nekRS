import numpy as np  # pylint: disable=import-error


def Locate(rank, nproc, datasize):
    extra = 0
    if rank == nproc - 1:
        extra = datasize % nproc
    num = datasize // nproc
    start = num * rank
    size = num + extra
    return start, size


class MPISetup(object):
    readargs = []
    size = 1
    rank = {"app": 0, "x": 0, "y": 0, "z": 0}

    def __init__(self, args, appID):

        self.nx = args.nx
        self.ny = args.ny
        self.nz = args.nz

        if not args.nompi:

            from mpi4py import MPI  # pylint: disable=import-error

            self.comm_app = MPI.COMM_WORLD.Split(appID, MPI.COMM_WORLD.Get_rank())
            self.size = self.comm_app.Get_size()
            self.rank["app"] = self.comm_app.Get_rank()
            if self.nx * self.ny * self.nz == 1:
                self.nx = self.size
            if self.size != (self.nx * self.ny * self.nz):
                raise ValueError("nx * ny * nz != num processes")

            if (self.ny > 1) and (self.nx > 1) and (self.nz > 1):
                comm_x = self.comm_app.Split(
                    self.rank["app"] % self.nx, self.rank["app"]
                )
            else:
                comm_x = self.comm_app.Split(
                    self.rank["app"] / self.nx, self.rank["app"]
                )
            comm_y = self.comm_app.Split(self.rank["app"] / self.ny, self.rank["app"])
            comm_z = self.comm_app.Split(self.rank["app"] / self.nz, self.rank["app"])

            self.rank["x"] = comm_x.Get_rank()
            self.rank["y"] = comm_y.Get_rank()
            self.rank["z"] = comm_z.Get_rank()

            self.readargs.append(self.comm_app)

        else:
            if self.nx != 1:
                raise ValueError("nx must = 1 without MPI")
            if self.ny != 1:
                raise ValueError("ny must = 1 without MPI")
            self.comm_app = None

    #            self.readargs.extend([args.xmlfile, "heat"])

    def Partition_3D_3D(self, fp, args):
        datashape = np.zeros(3, dtype=np.int64)
        start = np.zeros(3, dtype=np.int64)
        size = np.zeros(3, dtype=np.int64)

        # print("dir fp {0}".format(dir(fp)))
        var = fp.available_variables()
        #        print("var {0}".format(var))
        #        print('key')
        #        print('type {0}'.format(type(var)))
        # print("args varname {0}".format(args.varname))
        # print('var keys')
        #        print(var.keys())
        dshape = var[args.varname]["Shape"].split(",")
        for i in range(len(dshape)):
            datashape[i] = int(dshape[i])

        start[0], size[0] = Locate(self.rank["y"], self.ny, datashape[0])
        start[1], size[1] = Locate(self.rank["x"], self.nx, datashape[1])
        start[0], size[0] = Locate(self.rank["y"], self.ny, datashape[0])

        return start, size, datashape

    def Partition_3D_1D(self, fp, args):
        datashape = np.zeros(3, dtype=np.int64)
        start = np.zeros(3, dtype=np.int64)
        size = np.zeros(3, dtype=np.int64)

        # print("dir fp {0}".format(dir(fp)))
        var = fp.availablevariables()
        #        print("var {0}".format(var))
        #        print('key')
        #        print('type {0}'.format(type(var)))
        # print("args varname {0}".format(args.varname))
        # print('var keys')
        #        print(var.keys())
        dshape = var[args.varname]["Shape"].split(",")
        for i in range(len(dshape)):
            datashape[i] = int(dshape[i])

        start[0], size[0] = Locate(self.rank["app"], self.size, datashape[0])
        start[1], size[1] = (0, datashape[1])
        start[2], size[2] = (0, datashape[2])

        return start, size, datashape

    def Partition_2D_1D(self, fp, args):
        datashape = np.zeros(2, dtype=np.int64)
        start = np.zeros(2, dtype=np.int64)
        size = np.zeros(2, dtype=np.int64)

        var = fp.availablevariables()
        dshape = var[args.varname]["Shape"].split(",")
        for i in range(len(dshape)):
            datashape[i] = int(dshape[i])

        start[0], size[0] = Locate(self.rank["app"], self.size, datashape[0])
        start[1], size[1] = (0, datashape[1])

        return start, size, datashape
