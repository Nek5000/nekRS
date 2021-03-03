#!/usr/bin/env bash
set -e
for x in "$@"; do
  if [[ $x == *"="* ]]; then
    export ${x%%=*}="${x#*=}"
  fi
done

touch SIZE tmp.usr

export CFLAGS="${CFLAGS} -w -fPIC -mcmodel=medium"
export FFLAGS="${FFLAGS} -w -fPIC -fcray-pointer -mcmodel=medium -I../../"
${NEK5000_SOURCE_DIR}/bin/nekconfig -build-dep
${NEK5000_SOURCE_DIR}/bin/nekconfig

rm SIZE tmp.usr
