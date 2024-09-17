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

exec flake8 .
