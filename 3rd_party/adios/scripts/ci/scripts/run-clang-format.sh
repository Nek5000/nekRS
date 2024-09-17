#!/usr/bin/env bash

echo "---------- Begin ENV ----------"
env | sort
echo "----------  End ENV  ----------"

# If a source dir is given, then change to it.  Otherwise assume we're already
# in it
if [ -n "${SOURCE_DIR}" ]
then
  cd ${SOURCE_DIR}
fi

# Check C and C++ code with clang-format
find source plugins testing examples bindings  -regextype posix-extended -iregex '.*\.(h|c|cpp|tcc|cu)' ! -path "source/adios2/toolkit/derived/parser/pregen-source/*" | xargs clang-format -i
DIFF="$(git diff)"
if [ -n "${DIFF}" ]
then
  echo "clang-format:"
  echo "  Code format checks failed."
  echo "  Please run clang-format v16 your changes before committing:"
  echo "  You can use our CI image for this with: scripts/developer/run-clang-format.sh"
  echo "  The following changes are suggested:"
  echo "${DIFF}"
  echo "$(git diff --stat)"
  exit 1
fi


exit 0
