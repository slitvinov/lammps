#ifndef LAMMPS_MY_POOL_CHUNK_H
#define LAMMPS_MY_POOL_CHUNK_H 
namespace LAMMPS_NS {
template <class T> class MyPoolChunk {
 public:
  int ndatum;
  int nchunk;
  MyPoolChunk(int user_minchunk = 1, int user_maxchunk = 1, int user_nbin = 1,
              int user_chunkperpage = 1024, int user_pagedelta = 1);
  ~MyPoolChunk();
  T *get(int &index);
  T *get(int n, int &index);
  void put(int index);
  double size() const;
  int status() const { return errorflag; }
 private:
  int minchunk;
  int maxchunk;
  int nbin;
  int chunkperpage;
  int pagedelta;
  int binsize;
  int errorflag;
  T **pages;
  int *whichbin;
  int npage;
  int *freelist;
  int *freehead;
  int *chunksize;
  void allocate(int ibin);
};
}
#endif
