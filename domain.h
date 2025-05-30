#ifndef LMP_DOMAIN_H
#define LMP_DOMAIN_H
namespace LAMMPS_NS {
class Region;
class Domain : protected Pointers {
public:
  int box_exist;
  int dimension;
  int nonperiodic;
  int xperiodic, yperiodic, zperiodic;
  int periodicity[3];
  int boundary[3][2];
  int triclinic;
  double xprd, yprd, zprd;
  double xprd_half, yprd_half, zprd_half;
  double prd[3];
  double prd_half[3];
  double prd_lamda[3];
  double prd_half_lamda[3];
  double boxlo[3], boxhi[3];
  double boxlo_lamda[3], boxhi_lamda[3];
  double boxlo_bound[3], boxhi_bound[3];
  double corners[8][3];
  double minxlo, minxhi;
  double minylo, minyhi;
  double minzlo, minzhi;
  double sublo[3], subhi[3];
  double sublo_lamda[3], subhi_lamda[3];
  double xy, xz, yz;
  double h[6], h_inv[6];
  double h_rate[6], h_ratelo[3];
  int box_change;
  int box_change_size;
  int box_change_shape;
  int box_change_domain;
  int copymode;
  enum { NO_REMAP, X_REMAP, V_REMAP };
  Domain(class LAMMPS *);
  virtual void init();
  void set_initial_box(int expandflag = 1);
  virtual void set_global_box();
  virtual void set_local_box();
  virtual void reset_box();
  virtual void pbc();
  inline int minimum_image_check(double dx, double dy, double dz) {
    if (xperiodic && fabs(dx) > xprd_half)
      return 1;
    if (yperiodic && fabs(dy) > yprd_half)
      return 1;
    if (zperiodic && fabs(dz) > zprd_half)
      return 1;
    return 0;
  }

protected:
  double small[3];
};
} // namespace LAMMPS_NS
#endif
