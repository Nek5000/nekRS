#ifndef RANDOMVECTOR_HPP
#define RANDOMVECTOR_HPP

#include <mpi.h>
#include <random>
#include <vector>
#include <algorithm>
#include <iostream>
#include <vector>
#include <platform.hpp>

template <typename T> std::vector<T> randomVector(int N, T min = 0, T max = 1, bool deterministic = false)
{
  unsigned int seed;
  if(deterministic) {
    seed = platform->comm.mpiRank; 
  } else {
    // std::random_device is a non-deterministic uniform random bit generator,
    // although implementations are allowed to implement std::random_device using
    // a pseudo-random number engine if there is no support for non-deterministic 
    // random number generation.
    std::random_device rd;

    // poor man solution to ensure random seed values across all mpi ranks
    for(int i = 0; i < platform->comm.mpiRank; i++) {
      seed = rd(); 
    }
  }

  std::default_random_engine dev(seed);

  std::uniform_real_distribution<T> dist{min, max};
  auto gen = [&dist, &dev]() { return dist(dev); };

  std::vector<T> vec(N);
  std::generate(vec.begin(), vec.end(), gen);

  return vec;
}
#endif
