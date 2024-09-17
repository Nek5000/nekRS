/*
 * stream.cpp
 *
 *  Created on: Nov 2018
 *      Author: Norbert Podhorszki
 */

#include "adios2/common/ADIOSConfig.h"

#include "adiosStream.h"
#include "stream.h"

#ifdef ADIOS2_HAVE_HDF5_PARALLEL
#include "hdf5Stream.h"
#endif

ioGroup::~ioGroup() {}
Stream::Stream(const std::string &streamName, const adios2::Mode mode)
: streamName(streamName), mode(mode)
{
}
Stream::~Stream() {}

void Stream::fillArray(std::shared_ptr<VariableInfo> ov, double value)
{
    if (ov->type == "double")
    {
        double *a = reinterpret_cast<double *>(ov->data.data());
        for (size_t i = 0; i < ov->datasize / ov->elemsize; ++i)
        {
            a[i] = value;
        }
    }
    else if (ov->type == "float")
    {
        float v = static_cast<float>(value);
        float *a = reinterpret_cast<float *>(ov->data.data());
        for (size_t i = 0; i < ov->datasize / ov->elemsize; ++i)
        {
            a[i] = v;
        }
    }
    else if (ov->type == "int")
    {
        int v = static_cast<int>(value);
        int *a = reinterpret_cast<int *>(ov->data.data());
        for (size_t i = 0; i < ov->datasize / ov->elemsize; ++i)
        {
            a[i] = v;
        }
    }
}

std::shared_ptr<Stream> openStream(const std::string &streamName, std::shared_ptr<ioGroup> iogroup,
                                   const adios2::Mode mode, IOLib iolib, MPI_Comm comm,
                                   bool iotimer, size_t appid)
{
    std::shared_ptr<Stream> sp;
    switch (iolib)
    {
    case IOLib::ADIOS: {
        auto s = adiosStream(streamName, iogroup->adiosio, mode, comm, iotimer, appid);
        sp = std::make_shared<adiosStream>(s);
        break;
    }
#ifdef ADIOS2_HAVE_HDF5_PARALLEL
    case IOLib::HDF5: {
        auto s = hdf5Stream(streamName, mode, comm);
        sp = std::make_shared<hdf5Stream>(s);
        break;
    }
#endif
    }
    return sp;
}
