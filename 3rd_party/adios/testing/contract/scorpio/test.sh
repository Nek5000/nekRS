#!/bin/bash

set -x
set -e

# shellcheck disable=SC1091
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/setup.sh"

# Fail if is not set
build_dir="${build_dir:?}"
test_dir="${test_dir:?}"

mkdir -p "${test_dir}"
cd "${test_dir}"

mpiexec --oversubscribe -np 4 "${build_dir}/examples/adios/example3" -v

bpls -d example3_1.nc.bp.dir/example3_1.nc.bp.0 > 0.dump
diff -u 0.dump /opt/adios2/source/testing/contract/scorpio/0.dump

bpls -d example3_1.nc.bp.dir/example3_1.nc.bp.1 > 1.dump
diff -u 1.dump /opt/adios2/source/testing/contract/scorpio/1.dump

bpls -d example3_1.nc.bp.dir/example3_1.nc.bp.2 > 2.dump
diff -u 2.dump /opt/adios2/source/testing/contract/scorpio/2.dump

bpls -d example3_1.nc.bp.dir/example3_1.nc.bp.3 > 3.dump
diff -u 3.dump /opt/adios2/source/testing/contract/scorpio/3.dump
