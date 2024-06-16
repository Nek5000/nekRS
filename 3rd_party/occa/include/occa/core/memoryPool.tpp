#include <type_traits>

namespace occa {
  template <class T>
  occa::memory memoryPool::reserve(const dim_t entries) {
    if(std::is_void<T>::value) return reserve(entries, occa::dtype::byte); 
    return reserve(entries, occa::dtype::get<T>());
  }
}
