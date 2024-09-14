#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include "adios2/common/ADIOSConfig.h"
#include <atl.h>
#include <evpath.h>
#include <pthread.h>

#include "sst.h"

#include "cp_internal.h"
#include "ffs_marshal.h"
#include "zfp.h"

static zfp_type GetZFPType(int Type)
{
    if (Type == Int32)
    {
        return zfp_type_int32;
    }
    else if (Type == UInt32)
    {
        return zfp_type_int32;
    }
    else if (Type == Int64)
    {
        return zfp_type_int64;
    }
    else if (Type == Int64)
    {
        return zfp_type_int64;
    }
    else if (Type == UInt64)
    {
        return zfp_type_int64;
    }
    else if (Type == UInt64)
    {
        return zfp_type_int64;
    }
    else if (Type == Float)
    {
        return zfp_type_float;
    }
    else if (Type == Double)
    {
        return zfp_type_double;
    }
    else if (Type == Double)
    {
        return zfp_type_double;
    }
    return zfp_type_none;
}

extern int ZFPcompressionPossible(const int Type, int DimCount)
{
    return ((GetZFPType(Type) != zfp_type_none) && (DimCount < 4));
}

static zfp_field *GetZFPField(void *Data, const size_t DimCount, int Type, const size_t *Dimensions)
{
    zfp_type zfpType = GetZFPType(Type);
    zfp_field *field = NULL;

    if (zfpType == zfp_type_none)
        return NULL;

    switch (DimCount)
    {
    case 1:
        field = zfp_field_1d(Data, zfpType, Dimensions[0]);
        break;
    case 2:
        field = zfp_field_2d(Data, zfpType, Dimensions[0], Dimensions[1]);
        break;
    case 3:
        field = zfp_field_3d(Data, zfpType, Dimensions[0], Dimensions[1], Dimensions[2]);
        break;
    default:
        fprintf(stderr, "ZFP Compression not supported on %ld dimensional data\n", DimCount);
        exit(1);
    }
    return field;
}

static atom_t ZFPToleranceAtom = -1;
static atom_t ZFPRateAtom = -1;
static atom_t ZFPPrecisionAtom = -1;

zfp_stream *GetZFPStream(const size_t DimCount, int Type, attr_list Parameters)
{

    zfp_stream *zstream = zfp_stream_open(NULL);
    double Tolerance, Rate, Precision;

    if (ZFPToleranceAtom == -1)
    {
        ZFPToleranceAtom = attr_atom_from_string("ZFPTolernace");
        ZFPRateAtom = attr_atom_from_string("ZFPRate");
        ZFPPrecisionAtom = attr_atom_from_string("ZFPPrecision");
    }
    int hasTolerance = get_double_attr(Parameters, ZFPToleranceAtom, &Tolerance);
    int hasRate = get_double_attr(Parameters, ZFPRateAtom, &Rate);
    int hasPrecision = get_double_attr(Parameters, ZFPPrecisionAtom, &Precision);
    if (hasTolerance + hasRate + hasPrecision > 1)
        fprintf(stderr, "ERROR: zfp parameters Tolerance, "
                        "Rate, Precision are mutually "
                        "exclusive, only one of them is "
                        "mandatory, from "
                        "class CompressZfp Transform\n");

    if (hasTolerance)
    {
        zfp_stream_set_accuracy(zstream, Tolerance);
    }
    else if (hasRate)
    {
        zfp_stream_set_rate(zstream, Rate, GetZFPType(Type), DimCount, 0);
    }
    else if (hasPrecision)
    {
        zfp_stream_set_precision(zstream, Precision);
    }

    return zstream;
}

extern char *FFS_ZFPCompress(SstStream Stream, const size_t DimCount, int Type, void *Data,
                             const size_t *Count, size_t *ByteCountP)
{
    struct FFSWriterMarshalBase *Info = Stream->WriterMarshalData;
    void *bufferOut;
    zfp_field *field = GetZFPField(Data, DimCount, Type, Count);

    zfp_stream *stream = GetZFPStream(DimCount, Type, Info->ZFPParams);
    size_t maxSize = zfp_stream_maximum_size(stream, field);
    // associate bitstream
    bufferOut = malloc(maxSize);
    bitstream *bitstream = stream_open(bufferOut, maxSize);
    zfp_stream_set_bit_stream(stream, bitstream);
    zfp_stream_rewind(stream);
    size_t sizeOut = zfp_compress(stream, field);
    zfp_field_free(field);
    zfp_stream_close(stream);
    stream_close(bitstream);
    *ByteCountP = sizeOut;
    return bufferOut;
}

void *FFS_ZFPDecompress(SstStream Stream, const size_t DimCount, int Type, void *bufferIn,
                        const size_t sizeIn, const size_t *Dimensions, attr_list Parameters)
{
    zfp_field *in_field = GetZFPField(bufferIn, DimCount, Type, Dimensions);
    zfp_stream *stream = GetZFPStream(DimCount, Type, NULL);
    size_t maxSize = zfp_stream_maximum_size(stream, in_field);
    zfp_field_free(in_field);
    void *dataOut = malloc(maxSize);

    zfp_field *out_field = GetZFPField(dataOut, DimCount, Type, Dimensions);

    // associate bitstream
    bitstream *bitstream = stream_open(bufferIn, sizeIn);
    zfp_stream_set_bit_stream(stream, bitstream);
    zfp_stream_rewind(stream);

    int status = zfp_decompress(stream, out_field);

    if (!status)
        fprintf(stderr,
                "ERROR: zfp failed with status %d, in call to "
                "CompressZfp Decompress\n",
                status);

    zfp_field_free(out_field);
    zfp_stream_close(stream);
    stream_close(bitstream);

    return dataOut;
}
