#include "my_pool_chunk.h"
#include <cstdlib>
#if defined(LMP_INTEL) && !defined(LAMMPS_MEMALIGN) && !defined(_WIN32)
#define LAMMPS_MEMALIGN 64
#endif
using namespace LAMMPS_NS;
template <class T>
MyPoolChunk<T>::MyPoolChunk(int user_minchunk, int user_maxchunk, int user_nbin,
                            int user_chunkperpage, int user_pagedelta) {
  minchunk = user_minchunk;
  maxchunk = user_maxchunk;
  nbin = user_nbin;
  chunkperpage = user_chunkperpage;
  pagedelta = user_pagedelta;
  errorflag = 0;
  if (minchunk <= 0 || minchunk > maxchunk)
    errorflag = 1;
  if (user_nbin <= 0 || chunkperpage <= 0 || pagedelta <= 0)
    errorflag = 1;
  freehead = new int[nbin];
  chunksize = new int[nbin];
  if (!freehead || !chunksize)
    errorflag = 1;
  if (errorflag)
    return;
  binsize = (maxchunk - minchunk + 1) / nbin;
  if (minchunk + nbin * binsize <= maxchunk)
    binsize++;
  freelist = nullptr;
  for (int ibin = 0; ibin < nbin; ibin++) {
    freehead[ibin] = -1;
    chunksize[ibin] = minchunk + (ibin + 1) * binsize - 1;
    if (chunksize[ibin] > maxchunk)
      chunksize[ibin] = maxchunk;
  }
  ndatum = nchunk = 0;
  pages = nullptr;
  whichbin = nullptr;
  npage = 0;
}
template <class T> MyPoolChunk<T>::~MyPoolChunk() {
  delete[] freehead;
  delete[] chunksize;
  if (npage) {
    free(freelist);
    for (int i = 0; i < npage; i++)
      free(pages[i]);
    free(pages);
    free(whichbin);
  }
}
template <class T> T *MyPoolChunk<T>::get(int &index) {
  int ibin = nbin - 1;
  if (freehead[ibin] < 0) {
    allocate(ibin);
    if (errorflag) {
      index = -1;
      return nullptr;
    }
  }
  ndatum += maxchunk;
  nchunk++;
  index = freehead[ibin];
  int ipage = index / chunkperpage;
  int ientry = index % chunkperpage;
  freehead[ibin] = freelist[index];
  return &pages[ipage][ientry * chunksize[ibin]];
}
template <class T> T *MyPoolChunk<T>::get(int n, int &index) {
  if (n < minchunk || n > maxchunk) {
    errorflag = 3;
    index = -1;
    return nullptr;
  }
  int ibin = (n - minchunk) / binsize;
  if (freehead[ibin] < 0) {
    allocate(ibin);
    if (errorflag) {
      index = -1;
      return nullptr;
    }
  }
  ndatum += n;
  nchunk++;
  index = freehead[ibin];
  int ipage = index / chunkperpage;
  int ientry = index % chunkperpage;
  freehead[ibin] = freelist[index];
  return &pages[ipage][ientry * chunksize[ibin]];
}
template <class T> void MyPoolChunk<T>::put(int index) {
  if (index < 0)
    return;
  int ipage = index / chunkperpage;
  int ibin = whichbin[ipage];
  nchunk--;
  ndatum -= chunksize[ibin];
  freelist[index] = freehead[ibin];
  freehead[ibin] = index;
}
template <class T> void MyPoolChunk<T>::allocate(int ibin) {
  int oldpage = npage;
  npage += pagedelta;
  freelist = (int *)realloc(freelist, sizeof(int) * npage * chunkperpage);
  pages = (T **)realloc(pages, sizeof(T *) * npage);
  whichbin = (int *)realloc(whichbin, sizeof(int) * npage);
  if (!freelist || !pages) {
    errorflag = 2;
    return;
  }
  for (int i = oldpage; i < npage; i++) {
    whichbin[i] = ibin;
#if defined(LAMMPS_MEMALIGN)
    void *ptr;
    if (posix_memalign(&ptr, LAMMPS_MEMALIGN,
                       sizeof(T) * chunkperpage * chunksize[ibin]))
      errorflag = 2;
    pages[i] = (T *)ptr;
#else
    pages[i] = (T *)malloc(sizeof(T) * chunkperpage * chunksize[ibin]);
    if (!pages[i])
      errorflag = 2;
#endif
  }
  freehead[ibin] = oldpage * chunkperpage;
  for (int i = freehead[ibin]; i < npage * chunkperpage; i++)
    freelist[i] = i + 1;
  freelist[npage * chunkperpage - 1] = -1;
}
template <class T> double MyPoolChunk<T>::size() const {
  double bytes = (double)npage * chunkperpage * sizeof(int);
  bytes += (double)npage * sizeof(T *);
  bytes += (double)npage * sizeof(int);
  for (int i = 0; i < npage; ++i)
    bytes += (double)chunkperpage * chunksize[i] * sizeof(T);
  return bytes;
}
namespace LAMMPS_NS {
template class MyPoolChunk<int>;
template class MyPoolChunk<double>;
} // namespace LAMMPS_NS
