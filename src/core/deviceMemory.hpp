#ifndef deviceMemory_hpp_
#define deviceMemory_hpp_
#include "nekrsSys.hpp"

template <typename T> class DeviceMemoryAllocator
{
  using size_type = std::size_t;

public:
  static occa::memory malloc(size_type count)
  {
    if (std::is_same<T, std::byte>::value) {
      return platform->device.malloc(count);
    } else {
      return platform->device.malloc<T>(count);
    }
  }
};

template <typename T = std::byte, class Allocator = DeviceMemoryAllocator<T>> class deviceMemory
{
  using size_type = std::size_t;

public:
  deviceMemory() = default;

  explicit deviceMemory(size_type count) : occa_memory_{Allocator::malloc(count)} {}

  explicit deviceMemory(const std::vector<T> &values) : deviceMemory{values.size()}
  {
    copyFrom(values);
  }

  explicit deviceMemory(const occa::memory &occa_memory) : occa_memory_{occa_memory}
  {
    if (occa_memory_.byte_size() && occa_memory_.dtype() != occa::dtype::get<T>()) {
      throw std::runtime_error("data type does not match");
    }
  }

  explicit deviceMemory(occa::memory &&occa_memory) noexcept : occa_memory_{std::move(occa_memory)} {}

  deviceMemory(const deviceMemory &other) : occa_memory_{other.occa_memory_}
  {
    if (occa_memory_.byte_size() && occa_memory_.dtype() != occa::dtype::get<T>()) {
      throw std::runtime_error("data type does not match");
    }
  }

  deviceMemory(deviceMemory &&other) noexcept : occa_memory_{std::move(other.occa_memory_)} {}

  deviceMemory &operator=(const deviceMemory &rhs)
  {
    deviceMemory copy{rhs};
    swap(copy);
    return *this;
  }

  deviceMemory &operator=(deviceMemory &&rhs) noexcept
  {
    deviceMemory moved{std::move(rhs)};
    swap(moved);
    return *this;
  }

  ~deviceMemory() = default;

  operator const occa::memory &() const
  {
    return occa_memory_;
  }

  operator occa::memory &()
  {
    return occa_memory_;
  }

  operator occa::kernelArg() const
  {
    return occa_memory_;
  }

  void swap(deviceMemory &other) noexcept
  {
    occa::memory tmp = occa_memory_;
    occa_memory_ = other.occa_memory_;
    other.occa_memory_ = tmp;
  }

  size_type byte_size() const
  {
    return occa_memory_.byte_size();
  }

  size_type size() const
  {
    return occa_memory_.size();
  }

  size_type length() const
  {
    return occa_memory_.length();
  }

  bool isInitialized() const
  {
    return occa_memory_.isInitialized();
  }

  void clear()
  {
    occa_memory_.free();
  }

  void resize(size_type count)
  {
    if (count > size()) {
      if (size()) {
        clear();
      }
      occa_memory_ = Allocator::malloc(count);
    }
  }

  T *ptr()
  {
    return static_cast<T *>(occa_memory_.ptr());
  }

  const T *ptr() const
  {
    return static_cast<const T *>(occa_memory_.ptr());
  }

  deviceMemory slice(size_type offset, size_type count = 0) const
  {
    if (count) {
      return deviceMemory{occa_memory_.slice(offset, count)};
    } else {
      return deviceMemory{occa_memory_.slice(offset)};
    }
  }

  deviceMemory operator+(size_type offset) const
  {
    return slice(offset);
  };

  template <class A> void copyFrom(const deviceMemory<T, A> &src)
  {
    occa_memory_.copyFrom(src.occa_memory_);
  }

  template <class A>
  void copyFrom(const deviceMemory<T, A> &src,
                size_type count,
                size_type dest_offset = 0,
                size_type src_offset = 0)
  {
    occa_memory_.copyFrom(src.occa_memory_, count, dest_offset, src_offset);
  }

  void copyFrom(const std::vector<T> &src)
  {
    occa_memory_.copyFrom(src.data());
  }

  void copyFrom(const void *src, size_type count, size_type dest_offset = 0)
  {
    occa_memory_.copyFrom(src, count, dest_offset);
  }

  void copyFrom(const std::vector<T> &src, size_type count, size_type dest_offset = 0)
  {
    occa_memory_.copyFrom(src.data(), count, dest_offset);
  }

  void copyFrom(const occa::memory &src)
  {
    occa_memory_.copyFrom(src);
  }

  void copyFrom(const occa::memory &src, size_type count, size_type dest_offset = 0)
  {
    occa_memory_.copyFrom(src, count, dest_offset);
  }

  template <class A> void copyTo(deviceMemory<T, A> &dest) const
  {
    occa_memory_.copyTo(dest.occa_memory_);
  }

  template <class A>
  void
  copyTo(deviceMemory<T, A> &dest, size_type count, size_type dest_offset = 0, size_type src_offset = 0) const
  {
    occa_memory_.copyTo(dest.occa_memory_, count, dest_offset, src_offset);
  }

  void copyTo(std::vector<T> &dest) const
  {
    occa_memory_.copyTo(dest.data());
  }

  void copyTo(std::vector<T> &dest, size_type count, size_type src_offset = 0) const
  {
    occa_memory_.copyTo(dest.data(), count, src_offset);
  }

  void copyTo(occa::memory &dest) const
  {
    occa_memory_.copyTo(dest);
  }

  void copyTo(occa::memory &dest, size_type count, size_type src_offset = 0) const
  {
    occa_memory_.copyTo(dest, count, src_offset);
  }

private:
  occa::memory occa_memory_;
};

template <typename T> class DeviceMemoryPoolAllocator
{
  using size_type = std::size_t;

public:
  static occa::memory malloc(size_type count)
  {
    if (std::is_same<T, std::byte>::value) {
      return platform->deviceMemoryPool.reserve(count);
    } else {
      return platform->deviceMemoryPool.reserve<T>(count);
    }
  }
};

template <typename T = std::byte, class Allocator = DeviceMemoryPoolAllocator<T>>
class poolDeviceMemory : public deviceMemory<T>
{
  using size_type = std::size_t;

public:
  // Inherit base class constructors
  using deviceMemory<T>::deviceMemory;

  operator const deviceMemory<T> &() const
  {
    return deviceMemory<T>(occa_memory_);
  }

  operator deviceMemory<T> &()
  {
    return deviceMemory<T>(occa_memory_);
  }

private:
  occa::memory occa_memory_;
};

#endif
