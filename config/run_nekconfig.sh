#!/usr/bin/env bash
set -e
for x in "$@"; do
  if [[ $x == *"="* ]]; then
    export ${x%%=*}="${x#*=}"
  fi
done

export CFLAGS="${CFLAGS} -fPIC"
export FFLAGS="${FFLAGS} -fPIC -mcmodel=medium"
${NEK5000_SOURCE_DIR}/bin/nekconfig -build-dep
