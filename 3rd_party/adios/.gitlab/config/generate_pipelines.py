#!/usr/bin/env python3

# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#
# generate_pipeline.py
#
#   Author: Vicente Adolfo Bolea Sanchez <vicente.bolea@kitware.com>

from datetime import datetime
import argparse
import itertools
import requests
import time
import re
import urllib3
# Remove annoying warning about insecure connection (self-signed cert).
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


class skip_after_n_successes:
    def __init__(self, default_value, n):
        self.runs_max = n
        self.runs_current = 0
        self.default_value = default_value

    def __call__(self, fn, *args, **kwargs):
        if self.runs_current >= self.runs_max:
            return self.default_value

        ret = fn(*args, **kwargs)
        if ret:
            self.runs_current += 1
        return ret


def http_get_request(*args, **kwargs):
    kwargs['verify'] = False
    return requests.get(*args, **kwargs)


def request_as_list(url, *args, **kwargs):
    current_url = url
    body_json = []
    while current_url:
        response = http_get_request(current_url, *args, **kwargs)
        body_json += response.json()

        header = response.headers
        current_url = None
        if 'link' in header:
            links = re.search(
                r'(?<=\<)([\S]*)(?=>; rel="next")', header['link'], flags=re.IGNORECASE)
            if links is not None:
                current_url = links.group(0)

    return body_json


def add_timestamp(branch):
    date_str = branch['commit']['committed_date']
    # We ignore the TZ since Gitlab/GitHub always reports in UTC
    branch['dt'] = int(
        datetime.strptime(date_str.split(".")[0],
                          '%Y-%m-%dT%H:%M:%S').timestamp())
    return branch


def is_recent(branch):
    deadline_sec = int(time.time()) - (args.days * 86400)
    return branch['dt'] > deadline_sec


def has_no_status(branch):
    gh_commit_sha = branch['commit']['id']
    # Backported branches use the merge head
    if re.fullmatch(r'^pr\d+_.*$', branch['name']):
        gh_commit_sha = branch['commit']['parent_ids'][1]

    # Query GitHub for the status of this commit
    response = http_get_request(
        gh_url + '/commits/' + gh_commit_sha + '/status')
    if int(response.headers['x-ratelimit-remaining']) <= 0:
        raise ConnectionError(response.json())

    commit = response.json()
    if commit is None or 'sha' not in commit:
        return False

    for status in commit['statuses']:
        if status['context'] == args.gh_context:
            return False

    return True


parser = argparse.ArgumentParser(
    prog='generate_pipeline.py',
    description='Generate Dynamic pipelines for Gitlab')
parser.add_argument(
    '-u', '--gl-url', required=True,
    help='Base URL for Gitlab remote. Ex: https://code.olcf.ornl.gov/')
parser.add_argument(
    '-n', '--gh-name', required=True,
    help='Full name of the GitHub project. Ex: ornladios/ADIOS2')
parser.add_argument(
    '-p', '--project_id', required=True,
    help='Gitlab internal project ID of the project.')
parser.add_argument(
    '-c', '--gh-context', default='OLCF Crusher (Frontier)',
    help='Name of the status in GitHub (A.K.A context)')
parser.add_argument(
    '-d', '--days', type=int, default=1,
    help='How many days back to search for commits')
parser.add_argument(
    '-m', '--max', type=int, default=2,
    help='Maximum amount of pipelines computed')
parser.add_argument(
    '-f', '--template_file', required=True,
    help='Template file of the pipeline `{branch}` will be substituted')
args = parser.parse_args()


gl_url = args.gl_url + '/api/v4/projects/' + str(args.project_id)
gh_url = 'https://api.github.com/repos/' + args.gh_name

with open(args.template_file, 'r') as fd:
    template_str = fd.read()

    branches = request_as_list(gl_url + '/repository/branches')
    branches = [add_timestamp(branch) for branch in branches]
    branches = [b for b in branches if is_recent(b)]
    branches = sorted(branches, key=lambda x: x['dt'])

    # Skip running (and return true) has_no_status after returning True args.max times.
    # We need this not to hog the Github Rest API draconian ratelimit.
    run_n_times = skip_after_n_successes(default_value=False, n=args.max)
    branches = [b for b in branches if run_n_times(has_no_status, b)]

    for branch in branches:
        print(template_str.format(
            branch=branch['name'], commit=branch['commit']['id']))
