#ifndef LINSPACEVECTOR_HPP
#define LINSPACEVECTOR_HPP
#include <vector>
template <typename T> std::vector<T> linspace(T a, T b, int N)
{
  T dh = (b - a) / (N - 1);
  std::vector<T> points(N);

  T curr = a;
  for (auto &p : points) {
    p = curr;
    curr += dh;
  }
  return points;
}
#endif