#!/usr/bin/env bash
set -e
for x in "$@"; do
  if [[ $x == *"="* ]]; then
    export ${x%%=*}="${x#*=}"
  fi
done

make install
