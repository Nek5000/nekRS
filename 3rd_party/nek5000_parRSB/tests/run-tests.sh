#!/bin/bash

function info() {
  local GREEN='\033[0;92m'
  local NC='\033[0m'
  echo -e "${GREEN} $1 ${NC}"
}
export -f info

function error() {
  local RED='\033[0;31m'
  local NC='\033[0m'
  echo -e "${RED} $1 ${NC}"
}
export -f error

# Setup variables
tests=(`ls -I "*.okl" -I "*.sh" -I "*.log"`)

for t in "${tests[@]}"; do
  mpirun -use-hwthread-cpus -np 3 ./${t} >out.log 2>err.log
  wait $!
  if [ ! -s err.log ]; then
    info  "Test: ${t} ${b} ok."
  else
    error "Test: ${t} ${b} not ok."
  fi
done
