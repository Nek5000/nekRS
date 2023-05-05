#ifndef RANDOMVECTOR_HPP
#define RANDOMVECTOR_HPP

#include <mpi.h>
#include <random>
#include <vector>
#include <algorithm>
#include <iostream>
#include <vector>
template <typename T> std::vector<T> randomVector(int N)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::random_device rd;

  // poor man's solution seeding random_device
  unsigned int seed;
  for(int i=0; i<rank; i++) {
    seed = rd(); 
  }

  std::default_random_engine dev(seed);
  std::uniform_real_distribution<T> dist{0.0, 1.0};

  auto gen = [&dist, &dev]() { return dist(dev); };

  std::vector<T> vec(N);
  std::generate(vec.begin(), vec.end(), gen);

  return vec;
}
#endif
