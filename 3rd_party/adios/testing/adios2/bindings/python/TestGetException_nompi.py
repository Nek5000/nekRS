import numpy as np
import logging
import adios2.bindings as adios2

if __name__ == "__main__":
    __spec__ = None


def main():
    print("====================================")

    format = "%(asctime)s: %(message)s"
    logging.basicConfig(format=format, level=logging.INFO, datefmt="%H:%M:%S")
    writer()
    reader_fail()
    reader_success()


def writer():
    logging.info(" Writer: initiating writing")
    data = np.arange(3, dtype=np.float32)
    logging.info(f"        data dtypes: {data.dtype!s}")
    shape = data.shape
    count = shape
    start = (0,) * len(shape)
    logging.info(f"        data on writer side {data!s}")

    adios_io = adios2.ADIOS()
    IO = adios_io.DeclareIO("writer")

    writer = IO.Open("testdatafile", adios2.Mode.Write)
    writebuffer = IO.DefineVariable("np_data", data, shape, start, count, adios2.ConstantDims)
    if writebuffer:
        writer.BeginStep()
        writer.Put(writebuffer, data, adios2.Mode.Sync)
        writer.EndStep()
    else:
        raise ValueError("DefineVariable failed")

    writer.Close()

    logging.info(" Writer: writing finished")


def reader_fail():
    adios_io = adios2.ADIOS()
    io = adios_io.DeclareIO("reader")

    logging.info(" Reader: initiating erroneous reading ")
    reader = io.Open("testdatafile", adios2.Mode.Read)
    while True:
        stepStatus = reader.BeginStep()
        if stepStatus == adios2.StepStatus.OK:
            # inquire for variable
            recvar = io.InquireVariable("np_data")
            if recvar:
                got_exception = 0
                # determine the shape of the data that will be sent
                bufshape = recvar.Shape()
                typ = recvar.Type()
                logging.info("        Incoming variable type is " + typ)
                # allocate buffer for new numpy
                data = np.ones(bufshape)
                logging.info(f"        receive buffer  dtype: {data.dtype!s}")
                try:
                    reader.Get(recvar, data, adios2.Mode.Sync)
                except Exception as e:
                    got_exception = 1
                    logging.info(str(e))
                if not got_exception:
                    raise Exception("Didn't get an exception as expected")
                else:
                    logging.info("        Got expected exception")
            else:
                raise ValueError("InquireVariable failed")
        elif stepStatus == adios2.StepStatus.EndOfStream:
            break
        else:
            raise StopIteration(f"next step failed to initiate {stepStatus!s}")
        reader.EndStep()
    reader.Close()
    logging.info(
        " Reader: finished reading",
    )


def reader_success():
    adios_io = adios2.ADIOS()
    io = adios_io.DeclareIO("reader")

    logging.info(" Reader: initiating successful reading ")
    reader = io.Open("testdatafile", adios2.Mode.Read)
    while True:
        stepStatus = reader.BeginStep()
        if stepStatus == adios2.StepStatus.OK:
            # inquire for variable
            recvar = io.InquireVariable("np_data")
            if recvar:
                got_exception = 0
                # determine the shape of the data that will be sent
                bufshape = recvar.Shape()
                typ = recvar.Type()
                logging.info("        Incoming variable type is " + typ)
                # allocate buffer for new numpy
                data = np.ones(bufshape, dtype=np.float32)
                logging.info(f"        receive buffer  dtype: {data.dtype!s}")
                try:
                    reader.Get(recvar, data, adios2.Mode.Sync)
                except Exception as e:
                    got_exception = 1
                    logging.info(str(e))
                if got_exception:
                    raise Exception("        Got unexpected exception")
                else:
                    logging.info("        No exception as expected")
            else:
                raise ValueError("InquireVariable failed")
        elif stepStatus == adios2.StepStatus.EndOfStream:
            break
        else:
            raise StopIteration(f"next step failed to initiate {stepStatus!s}")
        reader.EndStep()
    reader.Close()
    logging.info(
        " Reader: finished reading",
    )


if __name__ == "__main__":
    main()
