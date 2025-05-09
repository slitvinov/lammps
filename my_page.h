#ifndef LAMMPS_MY_PAGE_H
#define LAMMPS_MY_PAGE_H
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
  T *vget() {
    if (index + maxchunk <= pagesize)
      return &page[index];
  }
  void vgot(int n) {
    ndatum += n;
    nchunk++;
    index += n;
  }
  void reset();
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
