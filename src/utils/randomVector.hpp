#include <random>
#include <vector>
#include <algorithm>
#include <iostream>
#include <vector>
template <typename T> std::vector<T> randomVector(int N)
{

  std::default_random_engine dev;
  std::uniform_real_distribution<T> dist{0.0, 1.0};

  auto gen = [&dist, &dev]() { return dist(dev); };

  std::vector<T> vec(N);
  std::generate(vec.begin(), vec.end(), gen);

  return vec;
}