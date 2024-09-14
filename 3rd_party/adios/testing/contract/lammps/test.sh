#!/bin/bash

set -x
set -e

# shellcheck disable=SC1091
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/setup.sh"

# Fail if is not set
install_dir="${install_dir:?}"
test_dir="${test_dir:?}"

mkdir -p "${test_dir}"
cd "${test_dir}"
cp -v /opt/adios2/source/testing/contract/lammps/{adios2_config.xml,check_results.sh,in.test} .


mpiexec -np 4 --oversubscribe "${install_dir}/bin/lmp" -in in.test

./check_results.sh
