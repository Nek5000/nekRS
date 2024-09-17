#!/bin/bash

set -ex

# Built the single el8 image
docker build --rm -f Dockerfile.ci-el8-intel -t adios2:ci-el8-intel .

# Tag images
docker tag adios2:ci-el8-intel ghcr.io/ornladios/adios2:ci-el8-oneapi
docker tag adios2:ci-el8-intel ghcr.io/ornladios/adios2:ci-el8-icc

# Push them
docker push ghcr.io/ornladios/adios2:ci-el8-oneapi
docker push ghcr.io/ornladios/adios2:ci-el8-icc
