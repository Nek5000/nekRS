#!/bin/bash

set -x
set -e

sudo /opt/spack/bin/spack install -v tau ~fortran ~papi ~pdt ~otf2
