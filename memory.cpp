#include "pointers.h"
#include "memory.h"
#if defined(LMP_INTEL) &&                                                      \
    ((defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)))
#ifndef LMP_INTEL_NO_TBB
#define LMP_USE_TBB_ALLOCATOR
#include "tbb/scalable_allocator.h"
#else
#include <cstring>
#include <malloc.h>
#endif
#endif
#if defined(LMP_INTEL) && !defined(LAMMPS_MEMALIGN) && !defined(_WIN32)
#define LAMMPS_MEMALIGN 64
#endif
using namespace LAMMPS_NS;
Memory::Memory(LAMMPS *lmp) : Pointers(lmp) {}
void *Memory::smalloc(bigint nbytes, const char *name) {
  if (nbytes == 0)
    return nullptr;
#if defined(LAMMPS_MEMALIGN)
  void *ptr;
  int retval = posix_memalign(&ptr, LAMMPS_MEMALIGN, nbytes);
  if (retval)
    ptr = nullptr;
#else
  void *ptr = malloc(nbytes);
#endif
  return ptr;
}
void *Memory::srealloc(void *ptr, bigint nbytes, const char *name) {
  if (nbytes == 0) {
    destroy(ptr);
    return nullptr;
  }
#if defined(LMP_INTEL_NO_TBB) && defined(LAMMPS_MEMALIGN) &&                   \
    (defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER))
  ptr = realloc(ptr, nbytes);
  uintptr_t offset = ((uintptr_t)(const void *)(ptr)) % LAMMPS_MEMALIGN;
  if (offset) {
    void *optr = ptr;
    ptr = smalloc(nbytes, name);
    memcpy(ptr, optr, MIN(nbytes, malloc_usable_size(optr)));
    free(optr);
  }
#else
  ptr = realloc(ptr, nbytes);
#endif
  return ptr;
}
void Memory::sfree(void *ptr) {
  if (ptr == nullptr)
    return;
  free(ptr);
}
void Memory::fail(const char *name) {
}
