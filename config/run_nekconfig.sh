#!/usr/bin/env bash
set -e
for x in "$@"; do
  if [[ $x == *"="* ]]; then
    export ${x%%=*}="${x#*=}"
  fi
done

${NEK5000_SOURCE_DIR}/bin/nekconfig -build-dep

touch SIZE tmp.usr
${NEK5000_SOURCE_DIR}/bin/nekconfig
mv makefile makefile.template
rm SIZE tmp.usr
