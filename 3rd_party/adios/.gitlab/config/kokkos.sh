#!/usr/bin/env bash
set -x

WORKDIR="$1"
VERSION="$2"

shift 2

if [ ! -d "$WORKDIR" ] || [ -z "$VERSION" ]
then
  echo "[E] missing args: Invoke as .gitlab/ci/config/kokkos.sh <WORKDIR> <VERSION> [extra_args]"
  exit 1
fi

# Build and install Kokkos
curl -L "https://github.com/kokkos/kokkos/archive/refs/tags/$VERSION.tar.gz" \
  | tar -C "$WORKDIR" -xzf -

cmake -S "$WORKDIR/kokkos-$VERSION" -B "$WORKDIR/kokkos_build" \
  "-DBUILD_SHARED_LIBS=ON" \
  "-DCMAKE_BUILD_TYPE:STRING=release" \
  "-DCMAKE_CXX_COMPILER_LAUNCHER=ccache" \
  "-DCMAKE_CXX_STANDARD:STRING=17" \
  "-DCMAKE_CXX_EXTENSIONS:BOOL=OFF" \
  "-DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON" \
  $*

cmake --build "$WORKDIR/kokkos_build"
cmake --install "$WORKDIR/kokkos_build"
