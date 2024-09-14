#!/usr/bin/env python3

import argparse
import glob
import sqlite3
import zlib
import yaml
from dataclasses import dataclass
from datetime import datetime
from os import chdir, getcwd, remove, stat
from os.path import exists, isdir, expanduser
from re import sub
from socket import getfqdn
from time import time_ns

# from adios2.adios2_campaign_manager import *

ADIOS_ACA_VERSION = "0.1"

@dataclass
class UserOption:
    adios_campaign_store: str = None
    hostname: str = None
    verbose: int = 0


def ReadUserConfig():
    path = expanduser("~/.config/adios2/adios2.yaml")
    opts = UserOption()
    try:
        doc = {}
        with open(path) as f:
            doc = yaml.safe_load(f)
        camp = doc.get("Campaign")
        if isinstance(camp, dict):
            for key, value in camp.items():
                if key == "campaignstorepath":
                    opts.adios_campaign_store = expanduser(value)
                if key == "hostname":
                    opts.hostname = value
                if key == "verbose":
                    opts.verbose = value
    except FileNotFoundError:
        None
    return opts


def SetupArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "command",
        help="Command: create/update/delete/info/list",
        choices=["create", "update", "delete", "info", "list"],
    )
    parser.add_argument(
        "campaign", help="Campaign name or path, with .aca or without", default=None, nargs="?"
    )
    parser.add_argument("--verbose", "-v", help="More verbosity", action="count", default=0)
    parser.add_argument(
        "--campaign_store", "-s", help="Path to local campaign store", default=None
    )
    parser.add_argument("--hostname", "-n", help="Host name unique for hosts in a campaign")
    parser.add_argument("-f", "--files", nargs="+", help="Add ADIOS files manually")
    args = parser.parse_args()

    # default values
    args.user_options = ReadUserConfig()

    if args.verbose == 0:
        args.verbose = args.user_options.verbose

    if args.campaign_store is None:
        args.campaign_store = args.user_options.adios_campaign_store

    if args.campaign_store is not None:
        while args.campaign_store[-1] == "/":
            args.campaign_store = args.campaign_store[:-1]

    if args.hostname is None:
        args.hostname = args.user_options.hostname

    args.CampaignFileName = args.campaign
    if args.campaign is not None:
        if not args.campaign.endswith(".aca"):
            args.CampaignFileName += ".aca"
        if (not exists(args.CampaignFileName) and
                not args.CampaignFileName.startswith("/") and
                args.campaign_store is not None):
            args.CampaignFileName = args.campaign_store + "/" + args.CampaignFileName

    if args.files is None:
        args.LocalCampaignDir = ".adios-campaign/"

    if args.verbose > 0:
        print(f"# Verbosity = {args.verbose}")
        print(f"# Command = {args.command}")
        print(f"# Campaign File Name = {args.CampaignFileName}")
        print(f"# Campaign Store = {args.campaign_store}")
    return args


def CheckCampaignStore(args):
    if args.campaign_store is not None and not isdir(args.campaign_store):
        print("ERROR: Campaign directory " + args.campaign_store + " does not exist", flush=True)
        exit(1)


def CheckLocalCampaignDir(args):
    if not isdir(args.LocalCampaignDir):
        print(
            "ERROR: Shot campaign data '" +
            args.LocalCampaignDir +
            "' does not exist. Run this command where the code was executed.",
            flush=True,
        )
        exit(1)


def IsADIOSDataset(dataset):
    if not isdir(dataset):
        return False
    if not exists(dataset + "/" + "md.idx"):
        return False
    if not exists(dataset + "/" + "data.0"):
        return False
    return True


def compressFile(f):
    compObj = zlib.compressobj()
    compressed = bytearray()
    blocksize = 1073741824  # 1GB #1024*1048576
    len_orig = 0
    len_compressed = 0
    block = f.read(blocksize)
    while block:
        len_orig += len(block)
        cBlock = compObj.compress(block)
        compressed += cBlock
        len_compressed += len(cBlock)
        block = f.read(blocksize)
    cBlock = compObj.flush()
    compressed += cBlock
    len_compressed += len(cBlock)

    return compressed, len_orig, len_compressed


def decompressBuffer(buf: bytearray):
    data = zlib.decompress(buf)
    return data


def AddFileToArchive(args: dict, filename: str, cur: sqlite3.Cursor, dsID: int):
    compressed = 1
    try:
        f = open(filename, "rb")
        compressed_data, len_orig, len_compressed = compressFile(f)

    except IOError:
        print(f"ERROR While reading file {filename}")
        return

    statres = stat(filename)
    ct = statres.st_ctime_ns

    cur.execute(
        "insert into bpfile "
        "(bpdatasetid, name, compression, lenorig, lencompressed, ctime, data) "
        "values (?, ?, ?, ?, ?, ?, ?) "
        "on conflict (bpdatasetid, name) do update "
        "set compression = ?, lenorig = ?, lencompressed = ?, ctime = ?, data = ?",
        (
            dsID,
            filename,
            compressed,
            len_orig,
            len_compressed,
            ct,
            compressed_data,
            compressed,
            len_orig,
            len_compressed,
            ct,
            compressed_data,
        ),
    )


def AddDatasetToArchive(hostID: int, dirID: int, dataset: str, cur: sqlite3.Cursor) -> int:
    statres = stat(dataset)
    ct = statres.st_ctime_ns
    select_cmd = (
        "select rowid from bpdataset "
        f"where hostid = {hostID} and dirid = {dirID} and name = '{dataset}'"
    )
    res = cur.execute(select_cmd)
    row = res.fetchone()
    if row is not None:
        rowID = row[0]
        print(
            f"Found dataset {dataset} in database on host {hostID} "
            f"in dir {dirID}, rowid = {rowID}"
        )
    else:
        print(f"Add dataset {dataset} to archive")
        curDS = cur.execute(
            "insert into bpdataset (hostid, dirid, name, ctime) values (?, ?, ?, ?)",
            (hostID, dirID, dataset, ct),
        )
        rowID = curDS.lastrowid
        # print(
        #     f"Inserted bpdataset {dataset} in database on host {hostID}"
        #     f" in dir {dirID}, rowid = {rowID}"
        # )
    return rowID


def ProcessFiles(args: dict, cur: sqlite3.Cursor, hostID: int, dirID: int):
    for entry in args.files:
        print(f"Process entry {entry}:")
        dsID = 0
        dataset = entry
        if IsADIOSDataset(dataset):
            dsID = AddDatasetToArchive(hostID, dirID, dataset, cur)
            cwd = getcwd()
            chdir(dataset)
            mdFileList = glob.glob("*md.*")
            profileList = glob.glob("profiling.json")
            files = mdFileList + profileList
            for f in files:
                AddFileToArchive(args, f, cur, dsID)
            chdir(cwd)
        else:
            print(f"WARNING: Dataset {dataset} is not an ADIOS dataset. Skip")


def GetHostName():
    host = getfqdn()
    if host.startswith("login"):
        host = sub("^login[0-9]*\\.", "", host)
    if host.startswith("batch"):
        host = sub("^batch[0-9]*\\.", "", host)
    if args.hostname is None:
        shorthost = host.split(".")[0]
    else:
        shorthost = args.user_options.hostname
    return host, shorthost


def AddHostName(longHostName, shortHostName):
    res = cur.execute('select rowid from host where hostname = "' + shortHostName + '"')
    row = res.fetchone()
    if row is not None:
        hostID = row[0]
        print(f"Found host {shortHostName} in database, rowid = {hostID}")
    else:
        curHost = cur.execute("insert into host values (?, ?)", (shortHostName, longHostName))
        hostID = curHost.lastrowid
        print(f"Inserted host {shortHostName} into database, rowid = {hostID}")
    return hostID


def MergeDBFiles(dbfiles: list):
    # read db files here
    result = list()
    for f1 in dbfiles:
        try:
            con = sqlite3.connect(f1)
        except sqlite3.Error as e:
            print(e)

        cur = con.cursor()
        try:
            cur.execute("select  * from bpfiles")
        except sqlite3.Error as e:
            print(e)
        record = cur.fetchall()
        for item in record:
            result.append(item[0])
        cur.close()
    return result


def AddDirectory(hostID: int, path: str) -> int:
    res = cur.execute(
        "select rowid from directory where hostid = " + str(hostID) + ' and name = "' + path + '"'
    )
    row = res.fetchone()
    if row is not None:
        dirID = row[0]
        print(f"Found directory {path} with hostID {hostID} in database, rowid = {dirID}")
    else:
        curDirectory = cur.execute("insert into directory values (?, ?)", (hostID, path))
        dirID = curDirectory.lastrowid
        print(f"Inserted directory {path} into database, rowid = {dirID}")
    return dirID


def Update(args: dict, cur: sqlite3.Cursor):
    longHostName, shortHostName = GetHostName()

    hostID = AddHostName(longHostName, shortHostName)

    rootdir = getcwd()
    dirID = AddDirectory(hostID, rootdir)
    con.commit()

    ProcessFiles(args, cur, hostID, dirID)

    con.commit()


def Create(args: dict, cur: sqlite3.Cursor):
    epoch = time_ns()
    cur.execute("create table info(id TEXT, name TEXT, version TEXT, ctime INT)")
    cur.execute(
        "insert into info values (?, ?, ?, ?)",
        ("ACA", "ADIOS Campaign Archive", ADIOS_ACA_VERSION, epoch),
    )
    cur.execute("create table host" + "(hostname TEXT PRIMARY KEY, longhostname TEXT)")
    cur.execute("create table directory" + "(hostid INT, name TEXT, PRIMARY KEY (hostid, name))")
    cur.execute(
        "create table bpdataset" +
        "(hostid INT, dirid INT, name TEXT, ctime INT" +
        ", PRIMARY KEY (hostid, dirid, name))"
    )
    cur.execute(
        "create table bpfile" +
        "(bpdatasetid INT, name TEXT, compression INT, lenorig INT" +
        ", lencompressed INT, ctime INT, data BLOB" +
        ", PRIMARY KEY (bpdatasetid, name))"
    )
    Update(args, cur)


def timestamp_to_datetime(timestamp: int) -> datetime:
    digits = len(str(int(timestamp)))
    t = float(timestamp)
    if digits > 18:
        t = t / 1000000000
    elif digits > 15:
        t = t / 1000000
    elif digits > 12:
        t = t / 1000
    return datetime.fromtimestamp(t)


def Info(args: dict, cur: sqlite3.Cursor):
    res = cur.execute("select id, name, version, ctime from info")
    info = res.fetchone()
    t = timestamp_to_datetime(info[3])
    print(f"{info[1]}, version {info[2]}, created on {t}")

    res = cur.execute("select rowid, hostname, longhostname from host")
    hosts = res.fetchall()
    for host in hosts:
        print(f"hostname = {host[1]}   longhostname = {host[2]}")
        res2 = cur.execute(
            'select rowid, name from directory where hostid = "' + str(host[0]) + '"'
        )
        dirs = res2.fetchall()
        for dir in dirs:
            print(f"    dir = {dir[1]}")
            res3 = cur.execute(
                'select rowid, name, ctime from bpdataset where hostid = "' +
                str(host[0]) +
                '" and dirid = "' +
                str(dir[0]) +
                '"'
            )
            bpdatasets = res3.fetchall()
            for bpdataset in bpdatasets:
                t = timestamp_to_datetime(bpdataset[2])
                print(f"        dataset = {bpdataset[1]}     created on {t}")


def List():
    path = args.campaign
    if path is None:
        if args.campaign_store is None:
            print("ERROR: Set --campaign_store for this command")
            return 1
        path = args.campaign_store
    else:
        while path[-1] == "/":
            path = path[:-1]

    # List the local campaign store
    acaList = glob.glob(path + "/**/*.aca", recursive=True)
    if len(acaList) == 0:
        print("There are no campaign archives in  " + path)
        return 2
    else:
        startCharPos = len(path) + 1
        for f in acaList:
            print(f[startCharPos:])
    return 0


def Delete():
    if exists(args.CampaignFileName):
        print(f"Delete archive {args.CampaignFileName}")
        remove(args.CampaignFileName)
        return 0
    else:
        print(f"ERROR: archive {args.CampaignFileName} does not exist")
        return 1


if __name__ == "__main__":
    args = SetupArgs()
    CheckCampaignStore(args)

    if args.command == "list":
        exit(List())

    if args.command == "delete":
        exit(Delete())

    if args.command == "create":
        print("Create archive")
        if exists(args.CampaignFileName):
            print(f"ERROR: archive {args.CampaignFileName} already exist")
            exit(1)
    elif args.command == "update" or args.command == "info":
        print(f"{args.command} archive")
        if not exists(args.CampaignFileName):
            print(f"ERROR: archive {args.CampaignFileName} does not exist")
            exit(1)

    con = sqlite3.connect(args.CampaignFileName)
    cur = con.cursor()

    if args.command == "info":
        Info(args, cur)
    else:
        if args.files is None:
            CheckLocalCampaignDir(args)
            # List the local campaign directory
            dbFileList = glob.glob(args.LocalCampaignDir + "/*.db")
            if len(dbFileList) == 0:
                print("There are no campaign data files in  " + args.LocalCampaignDir)
                exit(2)
            args.files = MergeDBFiles(dbFileList)

        if args.command == "create":
            Create(args, cur)
        elif args.command == "update":
            Update(args, cur)

    cur.close()
    con.close()
