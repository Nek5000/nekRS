/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * FileDaos.cpp file I/O using Daos I/O library
 *
 */
#include "FileDaos.h"
#include "adios2/helper/adiosLog.h"

#include <cstdio>      // remove
#include <cstdlib>     // malloc
#include <cstring>     // strerror
#include <errno.h>     // errno
#include <fcntl.h>     // open
#include <stddef.h>    // write output
#include <sys/stat.h>  // open, fstat
#include <sys/types.h> // open
#include <unistd.h>    // write, close

#include <algorithm>
#include <ios>
#include <iostream>

#include <chrono>

#include <daos.h>
#include <daos_fs.h>

// #include "adios2/helper/adiosFunctions.h"

#define CheckDAOSReturnCode(r) CheckDAOSReturnCodeF((r), __FILE__, __LINE__)
#define DefaultMaxDFSBatchSize 8589934592
namespace adios2
{
namespace transport
{

namespace
{
void GetUUIDFromEnv(const std::string &env, uuid_t &uuidValue)
{
    uuid_clear(uuidValue);

    if (const char *uuidEnv = std::getenv(env.c_str()))
    {
        if (!uuid_parse(uuidEnv, uuidValue))
        {
            helper::Throw<std::invalid_argument>(
                "Toolkit", "transport::file::FileDaos", "GetUUIDFromEnv",
                "Error: Unable to parse " + env + " environment variable");
        }
    }
    else
    {
        helper::Throw<std::invalid_argument>("Toolkit", "transport::file::FileDaos",
                                             "GetUUIDFromEnv",
                                             "Error: " + env + " environment variable not found");
    }
}

void CheckDAOSReturnCodeF(int rc, const char *file, int line)
{
    if (rc != 0)
    {
        std::string daosErr(d_errstr(rc));
        std::string theFile(file);
        helper::Throw<std::ios_base::failure>(
            "Toolkit", "transport::file::FileDaos", "CheckDAOSReturnCodeF",
            "ERROR: DAOS: " + daosErr + " " + theFile + " " + std::to_string(line));
    }
}
}

class FileDaos::Impl
{
public:
    Impl()
    {
        uuid_clear(UUID);
        try
        {
            GetUUIDFromEnv("ADIOS_DAOS_UUID", UUID);
        }
        catch (const std::invalid_argument &)
        {
        }

        uuid_clear(CUUID);
        try
        {
            GetUUIDFromEnv("ADIOS_DAOS_CUUID", CUUID);
        }
        catch (const std::invalid_argument &)
        {
        }

        if (const char *groupEnv = std::getenv("ADIOS_DAOS_GROUP"))
        {
            Group = groupEnv;
        }
    }

    ~Impl()
    {
        if (Obj)
        {
            int rc = dfs_release(Obj);
            CheckDAOSReturnCode(rc);
            Obj = nullptr;
        }
        if (Mount)
        {
            int rc = dfs_umount(Mount);
            CheckDAOSReturnCode(rc);
            Mount = nullptr;
        }
    }

    void Release()
    {
        if (Obj)
        {
            int rc = dfs_release(Obj);
            CheckDAOSReturnCode(rc);
            Obj = nullptr;
        }
        if (Mount)
        {
            int rc = dfs_umount(Mount);
            CheckDAOSReturnCode(rc);
            Mount = nullptr;
        }

        int rc;

        rc = daos_cont_close(coh, NULL);
        CheckDAOSReturnCode(rc);

        rc = daos_pool_disconnect(poh, NULL);
        CheckDAOSReturnCode(rc);
    }

    void InitMount(helper::Comm const &comm, const Mode openMode)
    {
        int rc;
        int poolFlags, contFlags, mountFlags;

        char uuid_c[100];
        char cuuid_c[100];
        uuid_unparse(UUID, uuid_c);
        uuid_unparse(CUUID, cuuid_c);
        std::string uuid_s{uuid_c};
        std::string cuuid_s{cuuid_c};

        if (Mount)
        {
            helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileDaos",
                                                  "InitMount", "Mount handle already exists");
        }

        rc = daos_init();
        CheckDAOSReturnCode(rc);

        switch (openMode)
        {
        case Mode::Write:
        case Mode::Append:
            poolFlags = DAOS_PC_RW;
            contFlags = DAOS_COO_RW;
            mountFlags = O_RDWR;
            break;
        case Mode::Read:
            poolFlags = DAOS_PC_RO;
            contFlags = DAOS_COO_RO;
            mountFlags = O_RDONLY;
            break;
        default:
            break;
        }

        if (singleProcess)
        {

#if DAOS_API_VERSION_MAJOR == 2
            rc = daos_pool_connect(uuid_c, Group.c_str(), poolFlags, &poh, NULL, NULL);
#else
            rc = daos_pool_connect(UUID, Group.c_str(), poolFlags, &poh, NULL, NULL);
#endif
            CheckDAOSReturnCode(rc);
#if DAOS_API_VERSION_MAJOR == 2
            rc = daos_cont_open(poh, cuuid_c, contFlags, &coh, NULL, NULL);
#else
            rc = daos_cont_open(poh, CUUID, contFlags, &coh, NULL, NULL);
#endif
            CheckDAOSReturnCode(rc);

            rc = dfs_mount(poh, coh, mountFlags, &Mount);
            CheckDAOSReturnCode(rc);
        }
        else
        {
            d_iov_t gHandles[3];
            std::vector<size_t> bufLen(3);

            memset(gHandles, 0, sizeof(d_iov_t) * 3);

            if (comm.Rank() == 0)
            {
#if DAOS_API_VERSION_MAJOR == 2
                rc = daos_pool_connect(uuid_c, Group.c_str(), poolFlags, &poh, NULL, NULL);
#else
                rc = daos_pool_connect(UUID, Group.c_str(), poolFlags, &poh, NULL, NULL);
#endif
                CheckDAOSReturnCode(rc);
#if DAOS_API_VERSION_MAJOR == 2
                rc = daos_cont_open(poh, cuuid_c, contFlags, &coh, NULL, NULL);
#else
                rc = daos_cont_open(poh, CUUID, contFlags, &coh, NULL, NULL);
#endif
                CheckDAOSReturnCode(rc);

                rc = dfs_mount(poh, coh, mountFlags, &Mount);
                CheckDAOSReturnCode(rc);

                rc = daos_pool_local2global(poh, &gHandles[0]);
                CheckDAOSReturnCode(rc);

                rc = daos_cont_local2global(coh, &gHandles[1]);
                CheckDAOSReturnCode(rc);

                rc = dfs_local2global(Mount, &gHandles[2]);
                CheckDAOSReturnCode(rc);

                if (gHandles[0].iov_buf == NULL)
                {
                    gHandles[0].iov_buf = malloc(gHandles[0].iov_buf_len);
                }
                if (gHandles[1].iov_buf == NULL)
                {
                    gHandles[1].iov_buf = malloc(gHandles[1].iov_buf_len);
                }
                if (gHandles[2].iov_buf == NULL)
                {
                    gHandles[2].iov_buf = malloc(gHandles[2].iov_buf_len);
                }
                rc = daos_pool_local2global(poh, &gHandles[0]);
                CheckDAOSReturnCode(rc);

                rc = daos_cont_local2global(coh, &gHandles[1]);
                CheckDAOSReturnCode(rc);

                rc = dfs_local2global(Mount, &gHandles[2]);
                CheckDAOSReturnCode(rc);

                for (size_t i = 0; i < 3; ++i)
                {
                    bufLen[i] = gHandles[i].iov_buf_len;
                }
            }
            comm.BroadcastVector(bufLen);
            size_t totalLen = 0;
            for (size_t i = 0; i < 3; ++i)
            {
                totalLen += bufLen[i];
                if (comm.Rank() != 0)
                {
                    gHandles[i].iov_buf_len = gHandles[i].iov_len = bufLen[i];
                }
            }

            std::vector<char> mergedBuf(totalLen);
            char *bufPtr = mergedBuf.data();

            if (comm.Rank() == 0)
            {
                for (size_t i = 0; i < 3; ++i)
                {
                    if (gHandles[i].iov_buf == NULL)
                    {
                        std::cout << "problem!" << std::endl;
                    }
                    std::memcpy(bufPtr, gHandles[i].iov_buf, gHandles[i].iov_buf_len);
                    bufPtr += gHandles[i].iov_buf_len;
                }
            }

            comm.BroadcastVector(mergedBuf);

            if (comm.Rank() != 0)
            {
                for (size_t i = 0; i < 3; ++i)
                {
                    gHandles[i].iov_buf = malloc(bufLen[i]);
                    std::memcpy(gHandles[i].iov_buf, bufPtr, bufLen[i]);
                    bufPtr += bufLen[i];
                }
                rc = daos_pool_global2local(gHandles[0], &poh);
                CheckDAOSReturnCode(rc);

                rc = daos_cont_global2local(poh, gHandles[1], &coh);
                CheckDAOSReturnCode(rc);

                rc = dfs_global2local(poh, coh, mountFlags, gHandles[2], &Mount);
                CheckDAOSReturnCode(rc);
            }

            for (size_t i = 0; i < 3; ++i)
            {
                free(gHandles[i].iov_buf);
                gHandles[i].iov_buf = NULL;
                gHandles[i].iov_buf_len = 0;
                gHandles[i].iov_len = 0;
            }
        }
    }

    bool singleProcess = false;
    uuid_t UUID;
    uuid_t CUUID;
    std::string Group;
    daos_handle_t poh;
    daos_handle_t coh;
    dfs_t *Mount = nullptr;
    dfs_obj_t *Obj = nullptr;
};

FileDaos::FileDaos(helper::Comm const &comm) : Transport("File", "Daos", comm), m_Impl(new Impl) {}

FileDaos::~FileDaos()
{
    // Cleanup will happen in the destructor of Impl
}

void FileDaos::SetParameters(const Params &params)
{
    // The parameters are first initialized by environment variabels if present
    // They are then overridden by config parameters if preset
    // If neither config mechanisms are available then an error is thrown

    {
        auto param = params.find("singleprocess");
        if (param != params.end())
        {
            std::string singleProc = param->second;
            if (singleProc == "true")
            {
                m_Impl->singleProcess = true;
            }
        }
    }
    // daos_pool_uuid
    {
        auto param = params.find("daos_pool_uuid");
        if (param != params.end())
        {
            if (!uuid_parse(param->second.c_str(), m_Impl->UUID))
            {
                helper::Throw<std::invalid_argument>("Toolkit", "transport::file::FileDaos",
                                                     "SetParameters",
                                                     "Unable to parse daos_pool_uuid parameter");
            }
        }
        if (uuid_is_null(m_Impl->UUID))
        {
            helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileDaos",
                                                  "SetParameters", "DAOS UUID is empty or not set");
        }
    }

    // daos_pool_group
    {
        auto param = params.find("daos_pool_group");
        if (param != params.end())
        {
            m_Impl->Group = param->second;
        }
        if (m_Impl->Group.empty())
        {
            helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileDaos",
                                                  "SetParameters",
                                                  "DAOS Group is empty or not set");
        }
    }

    // daos_cont_uuid
    {
        auto param = params.find("daos_cont_uuid");
        if (param != params.end())
        {
            if (!uuid_parse(param->second.c_str(), m_Impl->CUUID))
            {
                helper::Throw<std::invalid_argument>("Toolkit", "transport::file::FileDaos",
                                                     "SetParameters",
                                                     "Unable to parse daos_cont_uuid parameter");
            }
        }
        if (uuid_is_null(m_Impl->CUUID))
        {
            helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileDaos",
                                                  "SetParameters",
                                                  "DAOS CUUID is empty or not set");
        }
    }
}

void FileDaos::MkDir(const std::string &path)
{
    if (!m_Impl->Mount)
    {
        ProfilerStart("mount");
        m_Impl->InitMount(m_Comm, Mode::Write);
        ProfilerStop("mount");
    }
    struct stat stbuf;

    if (dfs_stat(m_Impl->Mount, NULL, path.c_str(), &stbuf) != 0)
    {
        if (!m_Impl->Mount)
        {
            std::cout << "m_Impl->Mount is NULL, problem!" << std::endl;
        }
        int rc = dfs_mkdir(m_Impl->Mount, NULL, path.c_str(), S_IWUSR | S_IRUSR, 0);
        CheckDAOSReturnCode(rc);
    }
    m_Impl->Release();
}

void FileDaos::WaitForOpen()
{
    if (m_IsOpening)
    {
        if (m_OpenFuture.valid())
        {
            m_DAOSOpenSucceed = m_OpenFuture.get();
        }
        m_IsOpening = false;
        CheckFile("couldn't open file " + m_Name + ", in call to Daos open");
        m_IsOpen = true;
    }
}

void FileDaos::Open(const std::string &name, const Mode openMode, const bool async,
                    const bool directio)
{

    int rc;
    m_Name = name;
    CheckName();

    if (!m_Impl->Mount)
    {
        m_Impl->InitMount(m_Comm, openMode);
    }
    m_OpenMode = openMode;

    dfs_obj_t *parent = NULL;
    std::string fileName;
    const auto lastPathSeparator(m_Name.find_last_of(PathSeparator));
    if (lastPathSeparator != std::string::npos)
    {
        const std::string dirName(m_Name.substr(0, lastPathSeparator));
        fileName = m_Name.substr(lastPathSeparator + 1, m_Name.length() - 1);
        rc = dfs_lookup_rel(m_Impl->Mount, NULL, dirName.c_str(), O_RDWR, &parent, NULL, NULL);
        CheckDAOSReturnCode(rc);
    }
    else
    {
        fileName = m_Name;
    }

    switch (m_OpenMode)
    {

    case Mode::Write:
        ProfilerStart("open");
        rc = dfs_open(/*DFS*/ m_Impl->Mount, /*PARENT*/ parent, fileName.c_str(), S_IFREG | S_IWUSR,
                      O_RDWR | O_CREAT, /*CID*/ 0,
                      /*chunksize*/ 0, NULL, &m_Impl->Obj);
        CheckDAOSReturnCode(rc);
        m_Errno = rc;
        ProfilerStop("open");
        break;

    case Mode::Append:
        ProfilerStart("open");
        rc = dfs_open(/*DFS*/ m_Impl->Mount, /*PARENT*/ parent, fileName.c_str(),
                      S_IFREG | S_IWUSR | S_IRUSR, O_RDWR | O_CREAT, /*CID*/ 0,
                      /*chunksize*/ 0, NULL, &m_Impl->Obj);
        CheckDAOSReturnCode(rc);
        m_Errno = rc;
        ProfilerStop("open");
        break;

    case Mode::Read:
        ProfilerStart("open");
        rc = dfs_open(/*DFS*/ m_Impl->Mount, /*PARENT*/ parent, fileName.c_str(), S_IFREG | S_IRUSR,
                      O_RDONLY, /*CID*/ 0,
                      /*chunksize*/ 0, NULL, &m_Impl->Obj);
        CheckDAOSReturnCode(rc);
        m_Errno = rc;
        ProfilerStop("open");
        break;

    default:
        CheckFile("unknown open mode for file " + m_Name + ", in call to Daos open");
    }
    rc = dfs_release(parent);
    CheckDAOSReturnCode(rc);
    if (rc == 0)
    {
        m_DAOSOpenSucceed = true;
    }
    else
    {
        m_DAOSOpenSucceed = false;
    }

    if (!m_IsOpening)
    {
        CheckFile("couldn't open file " + m_Name + ", in call to Daos open");
        m_IsOpen = true;
    }
}

void FileDaos::Write(const char *buffer, size_t size, size_t start)
{
    auto lf_Write = [&](const char *buffer, size_t size) {
        daos_size_t io_size = size;
        // daos_size_t written_size = 0;
        d_sg_list_t wsgl;
        d_iov_t iov;
        int rc = 0;

        d_iov_set(&iov, const_cast<char *>(buffer), size);
        wsgl.sg_nr = 1;
        wsgl.sg_nr_out = 1;
        wsgl.sg_iovs = &iov;
        wsgl.sg_iovs[0].iov_len = io_size;
        ProfilerStart("write");
        // errno = 0;

        // const auto writtenSize = write(m_FileDescriptor, buffer, size);
        // std::cout << "rank " << m_Comm.Rank() << ": start dfs_write..." <<
        // std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        rc = dfs_write(m_Impl->Mount, m_Impl->Obj, &wsgl, m_GlobalOffset, NULL);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
        std::cout << "rank " << m_Comm.Rank() << ": dfs_write took " << duration.count() << "s"
                  << std::endl;

        CheckDAOSReturnCode(rc);
        if (rc)
        {
            helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileDaos", "Write",
                                                  "couldn't write to file " + m_Name +
                                                      ", in call to Daos Write" + SysErrMsg());
        }
        // std::cout << "rank " << m_Comm.Rank() << ": dfs_write succeeded!" <<
        // std::endl;
        m_Errno = rc;
        ProfilerStop("write");
        m_GlobalOffset += size;
        //        while (written_size < size)
        //        {
        //            io_size = size - written_size;
        //            wsgl.sg_iovs[0].iov_len = io_size;
        //            ProfilerStart("write");
        //            // errno = 0;
        //
        //            // const auto writtenSize = write(m_FileDescriptor,
        //            buffer, size);
        //	    //std::cout << "rank " << m_Comm.Rank() << ": start
        // dfs_write..." << std::endl;
        //            rc = dfs_write(m_Impl->Mount, m_Impl->Obj, &wsgl,
        //            m_GlobalOffset,
        //                           NULL);
        //	    CheckDAOSReturnCode(rc);
        //            if (rc)
        //            {
        // helper::Throw<std::ios_base::failure>( "Toolkit",
        // "transport::file::FileDaos", "Write",
        //                    "couldn't write to file " + m_Name +
        //                    ", in call to Daos Write" + SysErrMsg());
        //            }
        //	    //std::cout << "rank " << m_Comm.Rank() << ": dfs_write
        // succeeded!" << std::endl;
        //            m_Errno = rc;
        //            ProfilerStop("write");
        //            /*if (writtenSize == -1)
        //            {
        //                if (errno == EINTR)
        //                {
        //                    continue;
        //                }
        // helper::Throw<std::ios_base::failure>( "Toolkit",
        // "transport::file::FileDaos", "Write",
        //                    "couldn't write to file " + m_Name +
        //                    ", in call to Daos Write" + SysErrMsg());
        //            }*/
        //
        //            buffer += io_size;
        //            written_size += io_size;
        //            m_GlobalOffset += io_size;
        //        }
    };

    // std::cout << "rank " << m_Comm.Rank() << ": start Write..." << std::endl;
    WaitForOpen();
    if (start != MaxSizeT)
    {
        m_GlobalOffset = start;
    }

    if (size > DefaultMaxDFSBatchSize)
    {
        // std::cout << "rank " << m_Comm.Rank() << ": size >
        // DefaultMaxDFSBatchSize" << std::endl;
        const size_t batches = size / DefaultMaxDFSBatchSize;
        const size_t remainder = size % DefaultMaxDFSBatchSize;

        size_t position = 0;
        for (size_t b = 0; b < batches; ++b)
        {
            lf_Write(&buffer[position], DefaultMaxDFSBatchSize);
            position += DefaultMaxDFSBatchSize;
        }
        lf_Write(&buffer[position], remainder);
    }
    else
    {
        // std::cout << "rank " << m_Comm.Rank() << ": size <=
        // DefaultMaxDFSBatchSize" << std::endl;
        lf_Write(buffer, size);
    }
}

void FileDaos::Read(char *buffer, size_t size, size_t start)
{
    auto lf_Read = [&](char *buffer, size_t size) {
        daos_size_t io_size = size;
        daos_size_t got_size;
        // daos_size_t read_size = 0;
        d_sg_list_t rsgl;
        d_iov_t iov;
        int rc = 0;

        d_iov_set(&iov, const_cast<char *>(buffer), size);
        rsgl.sg_nr = 1;
        rsgl.sg_nr_out = 1;
        rsgl.sg_iovs = &iov;
        rsgl.sg_iovs[0].iov_len = io_size;
        ProfilerStart("read");

        // std::cout << "rank " << m_Comm.Rank() << ": start dfs_read..." <<
        // std::endl;
        rc = dfs_read(m_Impl->Mount, m_Impl->Obj, &rsgl, m_GlobalOffset, &got_size, NULL);

        CheckDAOSReturnCode(rc);

        if (rc)
        {
            helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileDaos", "Read",
                                                  "couldn't read from file " + m_Name +
                                                      ", in call to Daos Read" + SysErrMsg());
        }
        // std::cout << "rank " << m_Comm.Rank() << ": dfs_read succeeded!" <<
        // std::endl;
        m_Errno = rc;
        ProfilerStop("read");
        m_GlobalOffset += size;
        //        while (read_size < size)
        //        {
        //            request_size = size - read_size;
        //            rsgl.sg_iovs[0].iov_len = request_size;
        //            ProfilerStart("read");
        //
        //	    //std::cout << "rank " << m_Comm.Rank() << ": start
        // dfs_read..." << std::endl;
        //            rc = dfs_read(m_Impl->Mount, m_Impl->Obj, &rsgl,
        //            m_GlobalOffset,
        //                          &got_size, NULL);
        //
        //	    CheckDAOSReturnCode(rc);
        //
        //            if (rc)
        //            {
        // helper::Throw<std::ios_base::failure>( "Toolkit",
        // "transport::file::FileDaos", "Read",
        //                "couldn't read
        //                from file " +
        //                                             m_Name + ", in call to
        //                                             Daos Read" +
        //                                             SysErrMsg());
        //            }
        //	    //std::cout << "rank " << m_Comm.Rank() << ": dfs_read
        // succeeded!" << std::endl;
        //            m_Errno = rc;
        //            ProfilerStop("read");
        //
        //            buffer += read_size;
        //            read_size += got_size;
        //            m_GlobalOffset += got_size;
        //        }
    };

    // std::cout << "rank " << m_Comm.Rank() << ": start Read..." << std::endl;
    WaitForOpen();

    if (start != MaxSizeT)
    {
        m_GlobalOffset = start;
    }

    if (size > DefaultMaxDFSBatchSize)
    {
        const size_t batches = size / DefaultMaxDFSBatchSize;
        const size_t remainder = size % DefaultMaxDFSBatchSize;

        size_t position = 0;
        for (size_t b = 0; b < batches; ++b)
        {
            lf_Read(&buffer[position], DefaultMaxDFSBatchSize);
            position += DefaultMaxDFSBatchSize;
        }
        lf_Read(&buffer[position], remainder);
    }
    else
    {
        lf_Read(buffer, size);
    }
}

size_t FileDaos::GetSize()
{
    WaitForOpen();
    daos_size_t file_size;
    int rc = 0;
    rc = dfs_get_size(m_Impl->Mount, m_Impl->Obj, &file_size);
    m_Errno = rc;
    return static_cast<size_t>(file_size);
}

void FileDaos::Flush() {}

void FileDaos::Close()
{
    WaitForOpen();
    ProfilerStart("close");
    int rc;
    rc = dfs_release(m_Impl->Obj);
    m_Impl->Obj = NULL;
    m_Errno = rc;
    ProfilerStop("close");

    if (rc)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileDaos", "Close",
                                              "couldn't close file " + m_Name +
                                                  ", in call to Daos IO close" + SysErrMsg());
    }

    m_IsOpen = false;
}

void FileDaos::Delete()
{
    WaitForOpen();
    if (m_IsOpen)
    {
        Close();
    }
    int rc;
    dfs_remove(m_Impl->Mount, NULL, m_Name.c_str(), true, NULL);
}

void FileDaos::CheckFile(const std::string hint) const
{
    if (!m_DAOSOpenSucceed)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileDaos", "CheckFile",
                                              "ERROR: " + hint + SysErrMsg());
    }
}

std::string FileDaos::SysErrMsg() const
{
    std::string desc(d_errdesc(m_Errno));
    return std::string(": errno = " + std::to_string(m_Errno) + ": " + desc);
}

void FileDaos::SeekToEnd() { m_GlobalOffset = GetSize(); }

void FileDaos::SeekToBegin() { m_GlobalOffset = 0; }

void FileDaos::Seek(const size_t start)
{
    if (start != MaxSizeT)
    {
        m_GlobalOffset = start;
    }
    else
    {
        SeekToEnd();
    }
}

void FileDaos::Truncate(const size_t length)
{
    helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileDaos", "Truncate",
                                          "Daos Truncate is not implemented yet");
}

} // end namespace transport
} // end namespace adios2