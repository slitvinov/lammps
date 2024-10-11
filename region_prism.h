#ifdef REGION_CLASS
RegionStyle(prism,RegPrism);
#else
#ifndef LMP_REGION_PRISM_H
#define LMP_REGION_PRISM_H 
#include "region.h"
namespace LAMMPS_NS {
class RegPrism : public Region {
  friend class CreateBox;
 public:
  RegPrism(class LAMMPS *, int, char **);
  ~RegPrism() override;
  int inside(double, double, double) override;
  int surface_interior(double *, double) override;
  int surface_exterior(double *, double) override;
 private:
  double xlo, xhi, ylo, yhi, zlo, zhi;
  double xy, xz, yz;
  double h[3][3], hinv[3][3];
  int dimension;
  double a[3], b[3], c[3];
  double clo[3], chi[3];
  double face[6][3];
  double corners[8][3];
  int tri[12][3];
  void find_nearest(double *, double &, double &, double &);
  int inside_tri(double *, double *, double *, double *, double *);
  double closest(double *, double *, double *, double);
};
}
#endif
#endif
