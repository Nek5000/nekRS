#!/usr/bin/env bash

set -e

CLEANUP_SSH_AGENT=0
CLEANUP_WORKDIR=0

function cleanup()
{
  echo "Cleanup:"
  if [ ${CLEANUP_WORKDIR} -eq 1 ]
  then
    echo "  Removing ${TMP_WORKDIR}"
    rm -rf "${TMP_WORKDIR}"
  fi
  if [ ${CLEANUP_SSH_AGENT} -eq 1 ]
  then
    echo "  Shutting down ssh-agent(${SSH_AGENT_PID})"
    eval $(ssh-agent -k)
  fi
}


if [ $# -ne 3 ]
then
  echo "Usage: $0 <github-repo> <gitlab-host> <gitlab-repo>"
  exit 1
fi

GITHUB_REPO=$1
GITLAB_HOST=$2
GITLAB_REPO=$3
GITLAB_URL="git@${GITLAB_HOST}:${GITLAB_REPO}.git"
BASEDIR="$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ -z "${GITLAB_SSH_KEY_BASE64}" ]
then
  echo "Error: GITLAB_SSH_KEY_BASE64 is empty"
  exit 1
fi

trap cleanup EXIT

# Start the ssh agent
eval $(ssh-agent -s)
CLEANUP_SSH_AGENT=1

# Add the necessary key
echo "${GITLAB_SSH_KEY_BASE64}" | base64 -d | tr -d '\r' | ssh-add -

git --version

export GIT_SSH_COMMAND="ssh -F /dev/null -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no"

# Setup the local repo with two remotes
TMP_WORKDIR=$(mktemp -d)
CLEANUP_WORKDIR=1
echo "Initializing local repo in ${TMP_WORKDIR}"
cd ${TMP_WORKDIR}
git init
git remote add github https://github.com/${GITHUB_REPO}.git
git remote add gitlab ${GITLAB_URL}

# Retrieve open PRs
OPEN_PR_BRANCHES="$(python3 ${BASEDIR}/github-list-prs.py -o ${GITHUB_REPO} | awk '{print "pr"$1"_"$2}' | sort -r)"
echo "Open PRs:"
for PR in ${OPEN_PR_BRANCHES}
do
  echo "  ${PR}"
done
echo

# Retrieve sync'd PRs 
SYNCD_PR_BRANCHES="$(echo $(git ls-remote gitlab github/pr* | awk '{print $2}' | sed 's|^refs/heads/github/||' | sort -r))"
echo "Syncd PRs:"
for PR in ${SYNCD_PR_BRANCHES}
do
  echo "  ${PR}"
done
echo

# Determine any closed PRs that are currently sync'd
SYNCD_CLOSED_PR_BRANCHES=""
if [ -n "${SYNCD_PR_BRANCHES}" ]
then
  for SPR in ${SYNCD_PR_BRANCHES}
  do
    if [ -n "${OPEN_PR_BRANCHES}" ]
    then
      IS_OPEN=0
      for OPR in ${OPEN_PR_BRANCHES}
      do
        if [ "${SPR}" = "${OPR}" ]
        then
          IS_OPEN=1
          break
        fi
      done
      if [ ${IS_OPEN} -eq 0 ]
      then
        SYNCD_CLOSED_PR_BRANCHES="${SYNCD_CLOSED_PR_BRANCHES} ${SPR}"
      fi
    fi
  done
fi
echo "Syncd Closed PRs:"
for PR in ${SYNCD_CLOSED_PR_BRANCHES}
do
  echo "  ${PR}"
done
echo

# Delete any sync'd PRs
if [ -n "${SYNCD_CLOSED_PR_BRANCHES}" ]
then
  CLOSED_REFSPECS=""
  for PR in ${SYNCD_CLOSED_PR_BRANCHES}
  do
    echo "Adding respec for closed ${PR}"
    CLOSED_REFSPECS="${CLOSED_REFSPECS} :github/${PR}"
  done
fi

# Sync open PRs to OLCF
if [ -n "${OPEN_PR_BRANCHES}" ]
then
  FETCH_REFSPECS=""
  OPEN_REFSPECS=""
  for PR in ${OPEN_PR_BRANCHES}
  do
    PR_NUM=$(expr "${PR}" : 'pr\([0-9]\+\)')
    FETCH_REFSPECS="${FETCH_REFSPECS} +refs/pull/${PR_NUM}/head:refs/remotes/github/${PR}"
    OPEN_REFSPECS="${OPEN_REFSPECS} github/${PR}:github/${PR}"
  done

  echo "Fetching GitHub refs for open PRs"
  git fetch -q github ${FETCH_REFSPECS}
  echo

  echo "Building local branches for open PRs"
  for PR in ${OPEN_PR_BRANCHES}
  do
    git branch -q github/${PR} github/${PR}
  done
  echo
fi

if [ -n "${CLOSED_REFSPECS}" -o -n "${OPEN_REFSPECS}" ]
then
  echo "Syncing PRs to GitLab"
  git push --porcelain -f gitlab ${CLOSED_REFSPECS} ${OPEN_REFSPECS}
  echo
fi
