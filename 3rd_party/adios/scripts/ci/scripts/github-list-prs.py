#!/usr/bin/env python3

import sys
import argparse
import urllib.request
import json

parser = argparse.ArgumentParser(description="List open PRs")
parser.add_argument("-t", "--token", help="GitHub access token")
group = parser.add_mutually_exclusive_group()
group.add_argument("-o", "--open", help="List open PRs",
                   action="store_true", default=False)
group.add_argument("-c", "--closed", help="List closed PRs",
                   action="store_true", default=False)
group.add_argument("-a", "--all", help="List all PRs",
                   action="store_true", default=False)
parser.add_argument("repo", help="GitHub repo (org/repo or user/repo)")
args = parser.parse_args()

if args.open:
    state = "open"
elif args.closed:
    state = "closed"
elif args.all:
    state = "all"
else:
    state = "open"

try:
    request = urllib.request.Request(
        'https://api.github.com/repos/%s/pulls?state=%s' % (args.repo, state))
    if args.token:
        request.add_header('Authorization', 'token %s' % args.token)
    response = urllib.request.urlopen(request)
except OSError:
    sys.exit(1)

prs = json.loads(response.read())

for pr in prs:
    print("%d %s" % (pr['number'], pr['head']['ref']))
