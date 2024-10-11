#include "memory.h"
#include "error.h"
#if defined(LMP_INTEL) && ((defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)))
#ifndef LMP_INTEL_NO_TBB
#define LMP_USE_TBB_ALLOCATOR 
#include "tbb/scalable_allocator.h"
#else
#include <cstring>
#if defined(__APPLE__)
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif
#endif
#endif
#if defined(LMP_INTEL) && !defined(LAMMPS_MEMALIGN) && !defined(_WIN32)
#define LAMMPS_MEMALIGN 64
#endif
using namespace LAMMPS_NS;
Memory::Memory(LAMMPS *lmp) : Pointers(lmp) {}
void *Memory::smalloc(bigint nbytes, const char *name)
{
  if (nbytes == 0) return nullptr;
#if defined(LAMMPS_MEMALIGN)
  void *ptr;
#if defined(LMP_USE_TBB_ALLOCATOR)
  ptr = scalable_aligned_malloc(nbytes, LAMMPS_MEMALIGN);
#else
  int retval = posix_memalign(&ptr, LAMMPS_MEMALIGN, nbytes);
  if (retval) ptr = nullptr;
#endif
#else
  void *ptr = malloc(nbytes);
#endif
  if (ptr == nullptr)
    error->one(FLERR,"Failed to allocate {} bytes for array {}",
                                 nbytes,name);
  return ptr;
}
void *Memory::srealloc(void *ptr, bigint nbytes, const char *name)
{
  if (nbytes == 0) {
    destroy(ptr);
    return nullptr;
  }
#if defined(LMP_USE_TBB_ALLOCATOR)
  ptr = scalable_aligned_realloc(ptr, nbytes, LAMMPS_MEMALIGN);
#elif defined(LMP_INTEL_NO_TBB) && defined(LAMMPS_MEMALIGN) && \
  (defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER))
  ptr = realloc(ptr, nbytes);
  uintptr_t offset = ((uintptr_t)(const void *)(ptr)) % LAMMPS_MEMALIGN;
  if (offset) {
    void *optr = ptr;
    ptr = smalloc(nbytes, name);
#if defined(__APPLE__)
    memcpy(ptr, optr, MIN(nbytes,malloc_size(optr)));
#elif defined(_WIN32) || defined(__MINGW32__)
    memcpy(ptr, optr, MIN(nbytes,_msize(optr)));
#else
    memcpy(ptr, optr, MIN(nbytes,malloc_usable_size(optr)));
#endif
    free(optr);
  }
#else
  ptr = realloc(ptr,nbytes);
#endif
  if (ptr == nullptr)
    error->one(FLERR,"Failed to reallocate {} bytes for array {}",
                                 nbytes,name);
  return ptr;
}
void Memory::sfree(void *ptr)
{
  if (ptr == nullptr) return;
  #if defined(LMP_USE_TBB_ALLOCATOR)
  scalable_aligned_free(ptr);
  #else
  free(ptr);
  #endif
}
void Memory::fail(const char *name)
{
  error->one(FLERR,"Cannot create/grow a vector/array of pointers for {}",name);
}
