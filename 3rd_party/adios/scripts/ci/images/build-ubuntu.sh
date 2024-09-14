#!/bin/bash

set -ex

# Build the base image
docker build --rm -f ./Dockerfile.ci-spack-ubuntu20.04-base -t ghcr.io/ornladios/adios2:ci-spack-ubuntu20.04-base .

# Which is also the gcc11 image
docker tag ghcr.io/ornladios/adios2:ci-spack-ubuntu20.04-base ghcr.io/ornladios/adios2:ci-spack-ubuntu20.04-gcc11

# Build the gcc8, gcc9, and gcc10 images
docker build --rm --build-arg GCC_VERSION=8 -f ./Dockerfile.ci-spack-ubuntu20.04-gcc -t ghcr.io/ornladios/adios2:ci-spack-ubuntu20.04-gcc8 .
docker build --rm --build-arg GCC_VERSION=9 -f ./Dockerfile.ci-spack-ubuntu20.04-gcc -t ghcr.io/ornladios/adios2:ci-spack-ubuntu20.04-gcc9 .
docker build --rm --build-arg GCC_VERSION=10 -f ./Dockerfile.ci-spack-ubuntu20.04-gcc -t ghcr.io/ornladios/adios2:ci-spack-ubuntu20.04-gcc10 .

# Build the clang6 and clang10 images
docker build --rm --build-arg CLANG_VERSION=6.0 -f ./Dockerfile.ci-spack-ubuntu20.04-clang -t ghcr.io/ornladios/adios2:ci-spack-ubuntu20.04-clang6 .
docker build --rm --build-arg CLANG_VERSION=10 -f ./Dockerfile.ci-spack-ubuntu20.04-clang -t ghcr.io/ornladios/adios2:ci-spack-ubuntu20.04-clang10 .

# Push images to github container registry
docker push ghcr.io/ornladios/adios2:ci-spack-ubuntu20.04-base
docker push ghcr.io/ornladios/adios2:ci-spack-ubuntu20.04-gcc8
docker push ghcr.io/ornladios/adios2:ci-spack-ubuntu20.04-gcc9
docker push ghcr.io/ornladios/adios2:ci-spack-ubuntu20.04-gcc10
docker push ghcr.io/ornladios/adios2:ci-spack-ubuntu20.04-gcc11
docker push ghcr.io/ornladios/adios2:ci-spack-ubuntu20.04-clang6
docker push ghcr.io/ornladios/adios2:ci-spack-ubuntu20.04-clang10
