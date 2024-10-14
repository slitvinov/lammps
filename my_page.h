#ifndef LAMMPS_MY_PAGE_H
#define LAMMPS_MY_PAGE_H
#include "lmptype.h"
namespace LAMMPS_NS {
struct HyperOneCoeff {
  double biascoeff;
  tagint tag;
};
template <class T> class MyPage {
public:
  int ndatum;
  int nchunk;
  MyPage();
  virtual ~MyPage();
  int init(int user_maxchunk = 1, int user_pagesize = 1024,
           int user_pagedelta = 1);
  T *get(int n = 1);
  T *vget() {
    if (index + maxchunk <= pagesize)
      return &page[index];
    ipage++;
    if (ipage == npage) {
      allocate();
      if (errorflag)
        return nullptr;
    }
    page = pages[ipage];
    index = 0;
    return &page[index];
  }
  void vgot(int n) {
    if (n > maxchunk)
      errorflag = 1;
    ndatum += n;
    nchunk++;
    index += n;
  }
  void reset();
  double size() const { return (double)npage * pagesize * sizeof(T); }
  int status() const { return errorflag; }

private:
  char padding[1024];
  T **pages;
  T *page;
  int npage;
  int ipage;
  int index;
  int maxchunk;
  int pagesize;
  int pagedelta;
  int errorflag;
  void allocate();
  void deallocate();
};
} // namespace LAMMPS_NS
#endif
