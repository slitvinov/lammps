#ifndef LMP_REGION_H
#define LMP_REGION_H
namespace LAMMPS_NS {
class Region : protected Pointers {
public:
  char *id, *style;
  Region **reglist;
  int interior;
  double xscale, yscale, zscale;
  double extent_xlo, extent_xhi;
  double extent_ylo, extent_yhi;
  double extent_zlo, extent_zhi;
  int bboxflag;
  int open_faces[6];
  double dx, dy, dz, theta;
  double v[3];
  double rpoint[3];
  double omega[3];
  double rprev;
  double xcenter[3];
  int vel_timestep;
  int nregion;
  Region(class LAMMPS *, int, char **);
  virtual void init();
  void prematch();
  int match(double, double, double);
  virtual int inside(double, double, double) = 0;

protected:
  void options(int, char **);
  double point[3], runit[3];

private:
  char *xstr, *ystr, *zstr, *tstr;
  int xvar, yvar, zvar, tvar;
  double axis[3];
};
} // namespace LAMMPS_NS
#endif
