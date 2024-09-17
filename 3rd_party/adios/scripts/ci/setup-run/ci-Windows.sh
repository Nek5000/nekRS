#!/bin/bash
set +u
# We do not use conda cmake, thus we have to hint cmake where is
# conda root dir.
export CMAKE_PREFIX_PATH="C:/Miniconda/Library;$CMAKE_PREFIX_PATH"
