#!/usr/bin/env python3
import json
import argparse
from os.path import exists


def SetupArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("--indent", "-i", type=int, default=2,
                        help="indent size, default=2")
    parser.add_argument("FILE", help="Name of the JSON file")
    args = parser.parse_args()

    # print(args)
    return args


def CheckFileName(args):
    if not exists(args.FILE):
        print("ERROR: File " + args.FILE + " does not exist", flush=True)
        exit(1)


if __name__ == "__main__":

    args = SetupArgs()
    CheckFileName(args)

    with open(args.FILE) as jsonfile:
        parsed = json.load(jsonfile)

    print(json.dumps(parsed, indent=args.indent, sort_keys=False))
