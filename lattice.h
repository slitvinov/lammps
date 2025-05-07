#ifndef LMP_LATTICE_H
#define LMP_LATTICE_H
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
};
} // namespace LAMMPS_NS
#endif
