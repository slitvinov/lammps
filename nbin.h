#ifndef LMP_NBIN_H
#define LMP_NBIN_H 
#include "pointers.h"
namespace LAMMPS_NS {
class NBin : protected Pointers {
 public:
  int istyle;
  bigint last_bin;
  double cutoff_custom;
  int nbinx, nbiny, nbinz;
  int mbins;
  int mbinx, mbiny, mbinz;
  int mbinxlo, mbinylo, mbinzlo;
  double binsizex, binsizey, binsizez;
  double bininvx, bininvy, bininvz;
  int *binhead;
  int *bins;
  int *atom2bin;
  int *nbinx_multi, *nbiny_multi, *nbinz_multi;
  int *mbins_multi;
  int *mbinx_multi, *mbiny_multi, *mbinz_multi;
  int *mbinxlo_multi, *mbinylo_multi, *mbinzlo_multi;
  double *binsizex_multi, *binsizey_multi, *binsizez_multi;
  double *bininvx_multi, *bininvy_multi, *bininvz_multi;
  int **binhead_multi;
  NBin(class LAMMPS *);
  ~NBin() override;
  void post_constructor(class NeighRequest *);
  virtual void copy_neighbor_info();
  virtual void bin_atoms_setup(int) = 0;
  virtual void setup_bins(int) = 0;
  virtual void bin_atoms() = 0;
  virtual double memory_usage() { return 0.0; }
 protected:
  int includegroup;
  double cutneighmin;
  double cutneighmax;
  int binsizeflag;
  double binsize_user;
  double *bboxlo, *bboxhi;
  int ncollections;
  double **cutcollectionsq;
  int dimension;
  int triclinic;
  int maxatom;
  int maxbin;
  int maxcollections;
  int *maxbins_multi;
  int coord2bin(double *);
  int coord2bin_multi(double *, int);
};
}
#endif
