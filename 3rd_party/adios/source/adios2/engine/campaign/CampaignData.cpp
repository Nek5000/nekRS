/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * CampaignData.cpp
 * Campaign data struct
 *
 *  Created on: May 16, 2023
 *      Author: Norbert Podhorszki pnorbert@ornl.gov
 */

#include "CampaignData.h"
#include "adios2/helper/adiosLog.h"
#include "adios2/helper/adiosString.h"

#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
#include <zlib.h>

namespace adios2
{
namespace core
{
namespace engine
{

/*
 * Data from processes to be recorded
 */

static int sqlcb_info(void *p, int argc, char **argv, char **azColName)
{
    CampaignData *cdp = reinterpret_cast<CampaignData *>(p);
    cdp->version = std::string(argv[0]);
    return 0;
};

static int sqlcb_host(void *p, int argc, char **argv, char **azColName)
{
    CampaignData *cdp = reinterpret_cast<CampaignData *>(p);
    CampaignHost ch;
    /*
    std::cout << "SQL: host record: ";
    for (int i = 0; i < argc; i++)
    {
        std::cout << azColName[i] << " = " << (argv[i] ? argv[i] : "NULL")
                  << std::endl;
    }
    std::cout << std::endl;
    */
    ch.hostname = std::string(argv[0]);
    ch.longhostname = std::string(argv[1]);
    cdp->hosts.push_back(ch);
    return 0;
};

static int sqlcb_directory(void *p, int argc, char **argv, char **azColName)
{
    CampaignData *cdp = reinterpret_cast<CampaignData *>(p);
    size_t hostid = helper::StringToSizeT(std::string(argv[0]), "SQL callback convert text to int");
    size_t hostidx = hostid - 1; // SQL rows start from 1, vector idx start from 0
    cdp->hosts[hostidx].directory.push_back(argv[1]);
    return 0;
};

static int sqlcb_bpdataset(void *p, int argc, char **argv, char **azColName)
{
    CampaignData *cdp = reinterpret_cast<CampaignData *>(p);
    CampaignBPDataset cds;
    size_t dsid = helper::StringToSizeT(std::string(argv[0]), "SQL callback convert text to int");
    size_t hostid = helper::StringToSizeT(std::string(argv[1]), "SQL callback convert text to int");
    size_t dirid = helper::StringToSizeT(std::string(argv[2]), "SQL callback convert text to int");
    cds.hostIdx = hostid - 1; // SQL rows start from 1, vector idx start from 0
    cds.dirIdx = dirid - 1;
    cds.name = argv[3];
    cdp->bpdatasets[dsid] = cds;
    return 0;
};

static int sqlcb_bpfile(void *p, int argc, char **argv, char **azColName)
{
    CampaignData *cdp = reinterpret_cast<CampaignData *>(p);
    CampaignBPFile cf;
    size_t dsid = helper::StringToSizeT(std::string(argv[0]), "SQL callback convert text to int");
    cf.bpDatasetIdx = dsid;
    cf.name = std::string(argv[1]);
    int comp = helper::StringTo<int>(std::string(argv[2]), "SQL callback convert text to int");
    cf.compressed = (bool)comp;
    cf.lengthOriginal =
        helper::StringToSizeT(std::string(argv[3]), "SQL callback convert text to int");
    cf.lengthCompressed =
        helper::StringToSizeT(std::string(argv[4]), "SQL callback convert text to int");
    cf.ctime = helper::StringTo<int64_t>(std::string(argv[5]), "SQL callback convert ctime to int");

    CampaignBPDataset &cds = cdp->bpdatasets[cf.bpDatasetIdx];
    cds.files.push_back(cf);
    return 0;
};

void ReadCampaignData(sqlite3 *db, CampaignData &cd)
{
    int rc;
    char *zErrMsg = 0;
    std::string sqlcmd;

    sqlcmd = "SELECT version FROM info";
    rc = sqlite3_exec(db, sqlcmd.c_str(), sqlcb_info, &cd, &zErrMsg);
    if (rc != SQLITE_OK)
    {
        std::cout << "SQL error: " << zErrMsg << std::endl;
        std::string m(zErrMsg);
        helper::Throw<std::invalid_argument>("Engine", "CampaignReader", "ReadCampaignData",
                                             "SQL error on reading info records:" + m);
        sqlite3_free(zErrMsg);
    }

    sqlcmd = "SELECT hostname, longhostname FROM host";
    rc = sqlite3_exec(db, sqlcmd.c_str(), sqlcb_host, &cd, &zErrMsg);
    if (rc != SQLITE_OK)
    {
        std::cout << "SQL error: " << zErrMsg << std::endl;
        std::string m(zErrMsg);
        helper::Throw<std::invalid_argument>("Engine", "CampaignReader", "ReadCampaignData",
                                             "SQL error on reading host records:" + m);
        sqlite3_free(zErrMsg);
    }

    sqlcmd = "SELECT hostid, name FROM directory";
    rc = sqlite3_exec(db, sqlcmd.c_str(), sqlcb_directory, &cd, &zErrMsg);
    if (rc != SQLITE_OK)
    {
        std::cout << "SQL error: " << zErrMsg << std::endl;
        std::string m(zErrMsg);
        helper::Throw<std::invalid_argument>("Engine", "CampaignReader", "ReadCampaignData",
                                             "SQL error on reading directory records:" + m);
        sqlite3_free(zErrMsg);
    }

    sqlcmd = "SELECT rowid, hostid, dirid, name FROM bpdataset";
    rc = sqlite3_exec(db, sqlcmd.c_str(), sqlcb_bpdataset, &cd, &zErrMsg);
    if (rc != SQLITE_OK)
    {
        std::cout << "SQL error: " << zErrMsg << std::endl;
        std::string m(zErrMsg);
        helper::Throw<std::invalid_argument>("Engine", "CampaignReader", "ReadCampaignData",
                                             "SQL error on reading bpdataset records:" + m);
        sqlite3_free(zErrMsg);
    }

    sqlcmd = "SELECT bpdatasetid, name, compression, lenorig, lencompressed, ctime "
             "FROM bpfile";
    rc = sqlite3_exec(db, sqlcmd.c_str(), sqlcb_bpfile, &cd, &zErrMsg);
    if (rc != SQLITE_OK)
    {
        std::cout << "SQL error: " << zErrMsg << std::endl;
        std::string m(zErrMsg);
        helper::Throw<std::invalid_argument>("Engine", "CampaignReader", "ReadCampaignData",
                                             "SQL error on reading bpfile records:" + m);
        sqlite3_free(zErrMsg);
    }
}

/* Decompress from in-memory source to file dest until stream ends or EOF.
   inf() returns Z_OK on success, Z_MEM_ERROR if memory could not be
   allocated for processing, Z_DATA_ERROR if the deflate data is
   invalid or incomplete, Z_VERSION_ERROR if the version of zlib.h and
   the version of the library linked do not match, or Z_ERRNO if there
   is an error reading or writing the files.
   http://www.zlib.net/zlib_how.html */
int inflateToFile(const unsigned char *source, const size_t blobsize, std::ofstream *dest)
{
    constexpr size_t CHUNK = 16777216;
    int ret;
    unsigned have;
    z_stream strm;

    std::vector<unsigned char> out(CHUNK);

    /* allocate inflate state */
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = 0;
    strm.next_in = Z_NULL;
    ret = inflateInit(&strm);
    if (ret != Z_OK)
        return ret;

    /* decompress until deflate stream ends or end of file */
    unsigned char *p = const_cast<unsigned char *>(source);
    uInt pos = 0;
    do
    {
        uInt CHUNK_SIZE = static_cast<uInt>(blobsize > CHUNK ? CHUNK : blobsize);
        strm.avail_in = CHUNK_SIZE;

        strm.next_in = p + pos;

        /* run inflate() on input until output buffer not full */
        do
        {
            strm.avail_out = CHUNK;
            strm.next_out = out.data();
            ret = inflate(&strm, Z_NO_FLUSH);
            switch (ret)
            {
            case Z_NEED_DICT:
                ret = Z_DATA_ERROR; /* and fall through */
            case Z_STREAM_ERROR:
            case Z_DATA_ERROR:
            case Z_MEM_ERROR:
                (void)inflateEnd(&strm);
                return ret;
            }
            have = CHUNK - strm.avail_out;
            dest->write(reinterpret_cast<char *>(out.data()), have);
            if (dest->bad())
            {
                helper::Throw<std::runtime_error>("Core", "Campaign", "Inflate",
                                                  "error writing file ");
            }

        } while (strm.avail_out == 0);
        pos += CHUNK_SIZE;
        /* done when inflate() says it's done */
    } while (ret != Z_STREAM_END);

    /* clean up and return */
    (void)inflateEnd(&strm);
    return ret == Z_STREAM_END ? Z_OK : Z_DATA_ERROR;
}

static int64_t timeToSec(int64_t ct)
{
    int64_t t;
    if (ct > 99999999999999999)
    {
        /* nanosec to sec */
        t = ct / 1000000000;
    }
    else if (ct > 99999999999999)
    {
        /* microsec to sec */
        t = ct / 1000000;
    }
    else if (ct > 99999999999)
    {
        /* millisec to sec */
        t = ct / 1000;
    }
    else
    {
        t = ct;
    }
    return t;
}

static bool isFileNewer(const std::string path, int64_t ctime)
{
    int result;
#ifdef _WIN32
    struct _stat s;
    result = _stat(path.c_str(), &s);
#else
    struct stat s;
    result = stat(path.c_str(), &s);
#endif
    if (result != 0)
    {
        return false;
    }

    int64_t ct = static_cast<int64_t>(s.st_ctime);
    int64_t ctSec = timeToSec(ct);
    int64_t ctimeSec = timeToSec(ctime);

    /*std::cout << "   Stat(" << path << "): size = " << s.st_size
              << " ct = " << ctSec << " ctime = " << ctimeSec << "\n";*/
    return (ctSec > ctimeSec);
}

void SaveToFile(sqlite3 *db, const std::string &path, const CampaignBPFile &bpfile)
{
    if (isFileNewer(path, bpfile.ctime))
    {
        return;
    }

    int rc;
    char *zErrMsg = 0;
    std::string sqlcmd;
    std::string id = std::to_string(bpfile.bpDatasetIdx);

    sqlite3_stmt *statement;
    sqlcmd =
        "SELECT data FROM bpfile WHERE bpdatasetid = " + id + " AND name = '" + bpfile.name + "'";
    // std::cout << "SQL statement: " << sqlcmd << "\n";
    rc = sqlite3_prepare_v2(db, sqlcmd.c_str(), static_cast<int>(sqlcmd.size()), &statement, NULL);
    if (rc != SQLITE_OK)
    {
        std::cout << "SQL error: " << zErrMsg << std::endl;
        std::string m(zErrMsg);
        helper::Throw<std::invalid_argument>("Engine", "CampaignReader", "SaveToFIle",
                                             "SQL error on reading info records:" + m);
        sqlite3_free(zErrMsg);
    }

    int result = 0;
    result = sqlite3_step(statement);
    if (result != SQLITE_ROW)
    {
        helper::Throw<std::invalid_argument>("Engine", "CampaignReader", "SaveToFIle",
                                             "Did not find record for :" + bpfile.name);
    }

    int iBlobsize = sqlite3_column_bytes(statement, 0);
    const void *p = sqlite3_column_blob(statement, 0);

    /*std::cout << "-- Retrieved from DB data of " << bpfile.name << " size = " << iBlobsize
              << " compressed = " << bpfile.compressed
              << " compressed size = " << bpfile.lengthCompressed
              << " original size = " << bpfile.lengthOriginal << " blob = " << p << "\n";*/

    size_t blobsize = static_cast<size_t>(iBlobsize);
    std::ofstream f;
    f.rdbuf()->pubsetbuf(0, 0);
    f.open(path, std::ios::out | std::ios::binary);
    if (bpfile.compressed)
    {
        const unsigned char *ptr = static_cast<const unsigned char *>(p);
        inflateToFile(ptr, blobsize, &f);
    }
    else
    {
        f.write(static_cast<const char *>(p), blobsize);
    }
    f.close();
}

} // end namespace engine
} // end namespace core
} // end namespace adios2
