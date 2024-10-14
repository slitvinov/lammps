#ifndef LMP_LATTICE_H
#define LMP_LATTICE_H
#include "pointers.h"
namespace LAMMPS_NS {
class Lattice : protected Pointers {
public:
  enum { NONE, SC, BCC, FCC, HCP, DIAMOND, SQ, SQ2, HEX, CUSTOM };
  int style;
  double xlattice, ylattice, zlattice;
  double a1[3], a2[3], a3[3];
  int nbasis;
  double **basis;
  Lattice(class LAMMPS *, int, char **);
  ~Lattice() override;
  void lattice2box(double &, double &, double &);
  void box2lattice(double &, double &, double &);
  void bbox(int, double, double, double, double &, double &, double &, double &,
            double &, double &);

private:
  double scale;
  double origin[3];
  int orientx[3];
  int orienty[3];
  int orientz[3];
  double primitive[3][3];
  double priminv[3][3];
  double rotaterow[3][3];
  double rotatecol[3][3];
  int orthogonal();
  int right_handed();
  int collinear();
  void setup_transform();
  void add_basis(double, double, double);
  double dot(double *, double *);
  void cross(double *, double *, double *);
};
} // namespace LAMMPS_NS
#endif
