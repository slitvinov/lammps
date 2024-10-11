#ifndef LMP_REGION_H
#define LMP_REGION_H 
#include "pointers.h"
namespace LAMMPS_NS {
class Region : protected Pointers {
 public:
  char *id, *style;
  Region **reglist;
  int interior;
  int scaleflag;
  double xscale, yscale, zscale;
  double extent_xlo, extent_xhi;
  double extent_ylo, extent_yhi;
  double extent_zlo, extent_zhi;
  int bboxflag;
  int varshape;
  int dynamic;
  int moveflag, rotateflag;
  int openflag;
  int open_faces[6];
  int copymode;
  struct Contact {
    double r;
    double delx, dely, delz;
    double radius;
    int iwall;
    int varflag;
  };
  Contact *contact;
  int cmax;
  int tmax;
  double dx, dy, dz, theta;
  double v[3];
  double rpoint[3];
  double omega[3];
  double rprev;
  double xcenter[3];
  double prev[5];
  int vel_timestep;
  int nregion;
  int size_restart;
  Region(class LAMMPS *, int, char **);
  ~Region() override;
  virtual void init();
  int dynamic_check();
  void prematch();
  int match(double, double, double);
  virtual void set_velocity();
  void velocity_contact(double *, double *, int);
  virtual int restart(char *, int &);
  virtual void length_restart_string(int &);
  virtual void reset_vel();
  virtual int inside(double, double, double) = 0;
  virtual int surface_interior(double *, double) = 0;
  virtual int surface_exterior(double *, double) = 0;
  virtual void shape_update() {}
  virtual void pretransform();
  virtual void set_velocity_shape() {}
  virtual void velocity_contact_shape(double *, double *) {}
 protected:
  void add_contact(int, double *, double, double, double);
  void options(int, char **);
  void point_on_line_segment(double *, double *, double *, double *);
  void forward_transform(double &, double &, double &);
  double point[3], runit[3];
 private:
  char *xstr, *ystr, *zstr, *tstr;
  int xvar, yvar, zvar, tvar;
  double axis[3];
  void inverse_transform(double &, double &, double &);
  void rotate(double &, double &, double &, double);
};
}
#endif
