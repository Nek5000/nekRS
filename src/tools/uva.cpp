#include <map>

#include <occa/core/base.hpp>
#include <occa/tools/misc.hpp>
#include <occa/tools/sys.hpp>
#include <occa/tools/uva.hpp>

namespace occa {
  ptrRangeMap uvaMap;
  memoryVector uvaStaleMemory;

  //---[ ptrRange ]---------------------
  ptrRange::ptrRange() :
    start(NULL),
    end(NULL) {}

  ptrRange::ptrRange(void *ptr, const udim_t bytes) :
    start((char*) ptr),
    end(((char*) ptr) + bytes) {}

  ptrRange::ptrRange(const ptrRange &other) :
    start(other.start),
    end(other.end) {}

  ptrRange& ptrRange::operator = (const ptrRange &other) {
    start = other.start;
    end   = other.end;

    return *this;
  }

  bool ptrRange::operator == (const ptrRange &other) const {
    return ((start < other.end) &&
            (end > other.start));
  }

  bool ptrRange::operator != (const ptrRange &other) const {
    return ((start >= other.end) ||
            (end <= other.start));
  }

  int operator < (const ptrRange &a, const ptrRange &b) {
    return ((a != b) && (a.start < b.start));
  }

  std::ostream& operator << (std::ostream& out,
                             const ptrRange &range) {
    out << '['
        << (void*) range.start << ", " << (void*) range.end
        << ']';
    return out;
  }
  //====================================


  //---[ UVA ]--------------------------
  occa::modeMemory_t* uvaToMemory(void *ptr) {
    if (!ptr) {
      return NULL;
    }
    ptrRangeMap::iterator it = uvaMap.find(ptr);
    return (it == uvaMap.end()) ? NULL : it->second;
  }

  bool isManaged(void *ptr) {
    occa::modeMemory_t *mem = uvaToMemory(ptr);
    if (mem) {
      return (mem->memInfo & uvaFlag::isManaged);
    }
    return false;
  }

  void startManaging(void *ptr) {
    occa::modeMemory_t *mem = uvaToMemory(ptr);
    if (mem) {
      mem->memInfo |= uvaFlag::isManaged;
    }
  }

  void stopManaging(void *ptr) {
    occa::modeMemory_t *mem = uvaToMemory(ptr);
    if (mem) {
      mem->memInfo &= ~uvaFlag::isManaged;
    }
  }

  void syncToDevice(void *ptr, const udim_t bytes) {
    occa::modeMemory_t *mem = uvaToMemory(ptr);
    if (mem) {
      syncMemToDevice(mem, bytes, ptrDiff(mem->uvaPtr, ptr));
    }
  }

  void syncToHost(void *ptr, const udim_t bytes) {
    occa::modeMemory_t *mem = uvaToMemory(ptr);
    if (mem) {
      syncMemToHost(mem, bytes, ptrDiff(mem->uvaPtr, ptr));
    }
  }

  void syncMemToDevice(occa::modeMemory_t *mem,
                       const udim_t bytes,
                       const udim_t offset) {

    if (mem) {
      occa::memory(mem).syncToDevice(bytes, offset);
    }
  }

  void syncMemToHost(occa::modeMemory_t *mem,
                     const udim_t bytes,
                     const udim_t offset) {

    if (mem) {
      occa::memory(mem).syncToHost(bytes, offset);
    }
  }

  bool needsSync(void *ptr) {
    occa::modeMemory_t *mem = uvaToMemory(ptr);
    return mem ? mem->isStale() : false;
  }

  void sync(void *ptr) {
    occa::modeMemory_t *mem = uvaToMemory(ptr);
    if (mem) {
      if (mem->inDevice()) {
        syncMemToHost(mem);
      } else {
        syncMemToDevice(mem);
      }
    }
  }

  void dontSync(void *ptr) {
    removeFromStaleMap(ptr);
  }

  void removeFromStaleMap(void *ptr) {
    ptrRangeMap::iterator it = uvaMap.find(ptr);
    if (it == uvaMap.end()) {
      return;
    }

    memory m(it->second);
    if (!m.uvaIsStale()) {
      return;
    }

    removeFromStaleMap(m.getModeMemory());
  }

  void removeFromStaleMap(modeMemory_t *mem) {
    if (!mem) {
      return;
    }

    occa::memory m(mem);
    const size_t staleEntries = uvaStaleMemory.size();

    for (size_t i = 0; i < staleEntries; ++i) {
      if (uvaStaleMemory[i] == mem) {
        m.uvaMarkFresh();
        uvaStaleMemory.erase(uvaStaleMemory.begin() + i);
        break;
      }
    }
  }

  void freeUvaPtr(void *ptr) {
    if (!ptr) {
      return;
    }
    modeMemory_t *modeMemory = uvaToMemory(ptr);
    if (modeMemory) {
      occa::memory(modeMemory).free();
      return;
    }
    OCCA_FORCE_ERROR("Freeing a non-uva pointer");
  }
  //====================================
}
