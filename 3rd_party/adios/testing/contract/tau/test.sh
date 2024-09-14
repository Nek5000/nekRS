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

TAU=$(spack location -i tau)/bin/tau_exec

mpiexec -np 2 "${TAU}" "${install_dir}/bin/adios2_basics_variablesShapes"

[ ! -f profile.0.0.0 ] || [ ! -s profile.0.0.0 ] && { echo "Error: file profile.0.0.0 not found or empty"; exit 1; }
[ ! -f profile.1.0.0 ] || [ ! -s profile.1.0.0 ] && { echo "Error: file profile.1.0.0 not found or empty"; exit 1; }

cat profile.0.0.0
cat profile.1.0.0
