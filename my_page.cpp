#include "lmptype.h"
#include "my_page.h"
#if defined(LMP_INTEL) && !defined(LAMMPS_MEMALIGN) && !defined(_WIN32)
#define LAMMPS_MEMALIGN 64
#endif
using namespace LAMMPS_NS;
template <class T>
MyPage<T>::MyPage()
    : ndatum(0), nchunk(0), pages(nullptr), page(nullptr), npage(0), ipage(-1),
      index(-1), maxchunk(-1), pagesize(-1), pagedelta(1), errorflag(0){};
template <class T> MyPage<T>::~MyPage() { deallocate(); }
template <class T>
int MyPage<T>::init(int user_maxchunk, int user_pagesize, int user_pagedelta) {
  maxchunk = user_maxchunk;
  pagesize = user_pagesize;
  pagedelta = user_pagedelta;
  if (maxchunk <= 0 || pagesize <= 0 || pagedelta <= 0)
    return 1;
  if (maxchunk > pagesize)
    return 1;
  deallocate();
  allocate();
  if (errorflag)
    return 2;
  reset();
  return 0;
}
template <class T> void MyPage<T>::reset() {
  ndatum = nchunk = 0;
  index = ipage = 0;
  page = (pages != nullptr) ? pages[ipage] : nullptr;
  errorflag = 0;
}
template <class T> void MyPage<T>::allocate() {
  npage += pagedelta;
  pages = (T **)realloc(pages, npage * sizeof(T *));
  if (!pages) {
    errorflag = 2;
    return;
  }
  for (int i = npage - pagedelta; i < npage; i++) {
#if defined(LAMMPS_MEMALIGN)
    void *ptr;
    if (posix_memalign(&ptr, LAMMPS_MEMALIGN, pagesize * sizeof(T)))
      errorflag = 2;
    pages[i] = (T *)ptr;
#else
    pages[i] = (T *)malloc(pagesize * sizeof(T));
    if (!pages[i])
      errorflag = 2;
#endif
  }
}
template <class T> void MyPage<T>::deallocate() {
  reset();
  for (int i = 0; i < npage; i++)
    free(pages[i]);
  free(pages);
  pages = nullptr;
  npage = 0;
}
namespace LAMMPS_NS {
template class MyPage<int>;
template class MyPage<long>;
template class MyPage<long long>;
template class MyPage<double>;
template class MyPage<HyperOneCoeff>;
} // namespace LAMMPS_NS
