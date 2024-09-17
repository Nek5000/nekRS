#!/bin/sh

if command -v apt-get >/dev/null
then
  apt-get update
  apt-get install -y sudo
  apt-get clean
elif command -v yum >/dev/null
then
  yum install -y sudo
  yum clean all
fi 

useradd -m -s /bin/bash adios
echo "adios ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/adios
