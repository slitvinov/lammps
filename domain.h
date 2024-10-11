#ifndef LMP_DOMAIN_H
#define LMP_DOMAIN_H 
#include "pointers.h"
#include <cmath>
#include <map>
#include <unordered_set>
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
  int deform_flag;
  int deform_vremap;
  int deform_groupbit;
  class Lattice *lattice;
  int copymode;
  enum { NO_REMAP, X_REMAP, V_REMAP };
  typedef Region *(*RegionCreator)(LAMMPS *, int, char **);
  typedef std::map<std::string, RegionCreator> RegionCreatorMap;
  RegionCreatorMap *region_map;
  Domain(class LAMMPS *);
  ~Domain() override;
  virtual void init();
  void set_initial_box(int expandflag = 1);
  virtual void set_global_box();
  virtual void set_lamda_box();
  virtual void set_local_box();
  virtual void reset_box();
  virtual void pbc();
  void subbox_too_small_check(double);
  void minimum_image(double &, double &, double &) const;
  void minimum_image(double *delta) const { minimum_image(delta[0], delta[1], delta[2]); }
  int closest_image(int, int);
  int closest_image(const double *const, int);
  void closest_image(const double *const, const double *const, double *const);
  void remap(double *, imageint &);
  void remap(double *);
  void remap_near(double *, double *);
  void unmap(double *, imageint);
  void unmap(const double *, imageint, double *);
  void image_flip(int, int, int);
  int ownatom(int, double *, imageint *, int);
  void set_lattice(int, char **);
  void add_region(int, char **);
  void delete_region(Region *);
  void delete_region(const std::string &);
  Region *get_region_by_id(const std::string &) const;
  const std::vector<Region *> get_region_by_style(const std::string &) const;
  const std::vector<Region *> get_region_list();
  void set_boundary(int, char **, int);
  void print_box(const std::string &);
  void boundary_string(char *);
  virtual void lamda2x(int);
  virtual void x2lamda(int);
  virtual void lamda2x(double *, double *);
  virtual void x2lamda(double *, double *);
  int inside(double *);
  int inside_nonperiodic(double *);
  void x2lamda(double *, double *, double *, double *);
  void bbox(double *, double *, double *, double *);
  void box_corners();
  void subbox_corners();
  void lamda_box_corners(double *, double *);
  inline int minimum_image_check(double dx, double dy, double dz)
  {
    if (xperiodic && fabs(dx) > xprd_half) return 1;
    if (yperiodic && fabs(dy) > yprd_half) return 1;
    if (zperiodic && fabs(dz) > zprd_half) return 1;
    return 0;
  }
 protected:
  double small[3];
  std::unordered_set<Region *> regions;
};
}
#endif
