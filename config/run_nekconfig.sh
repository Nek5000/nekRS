#!/usr/bin/env bash
set -xe
for x in "$@"; do
  if [[ $x == *"="* ]]; then
    export ${x%%=*}="${x#*=}"
  fi
done

export NEK_SOURCE_ROOT=$NEK5000_SOURCE_DIR

${NEK5000_SOURCE_DIR}/bin/nekconfig -build-dep

touch SIZE tmp.usr
${NEK5000_SOURCE_DIR}/bin/nekconfig
mv makefile makefile.template
rm SIZE tmp.usr
