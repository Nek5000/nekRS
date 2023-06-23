#include "device.hpp"

template <class T>
occa::memory device_t::mallocHost(size_t entries)
{
  occa::properties props;
  props["host"] = true;

  void *buffer = std::calloc(entries, sizeof(T));
  occa::memory h_scratch = _device.malloc<T>(entries, buffer, props);
  std::free(buffer);
  return h_scratch;
}

template <class T>
occa::memory device_t::malloc(size_t entries, const occa::memory& src)
{
  return _device.malloc<T>(entries, src);
}

template <class T>
occa::memory device_t::malloc(size_t entries, const void *src)
{
  return _device.malloc<T>(entries, src);
}

template <class T>
occa::memory device_t::malloc(size_t entries)
{
  void *buffer = std::calloc(entries, sizeof(T));
  occa::memory o_returnValue = _device.malloc<T>(entries, buffer);
  std::free(buffer);
  return o_returnValue;
}
