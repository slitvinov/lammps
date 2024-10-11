#ifndef LMP_GROUP_H
#define LMP_GROUP_H 
#include "pointers.h"
#include <map>
namespace LAMMPS_NS {
class Region;
class Group : protected Pointers {
 public:
  int ngroup;
  char **names;
  int *bitmask;
  int *inversemask;
  int *dynamic;
  Group(class LAMMPS *);
  ~Group() override;
  void assign(int, char **);
  void assign(const std::string &);
  void create(const std::string &, int *);
  int find(const std::string &);
  int find_or_create(const char *);
  bigint count_all();
  bigint count(int);
  bigint count(int, Region *);
  double mass(int);
  double mass(int, Region *);
  double charge(int);
  double charge(int, Region *);
  void bounds(int, double *);
  void bounds(int, double *, Region *);
  void xcm(int, double, double *);
  void xcm(int, double, double *, Region *);
  void vcm(int, double, double *);
  void vcm(int, double, double *, Region *);
  void fcm(int, double *);
  void fcm(int, double *, Region *);
  double ke(int);
  double ke(int, Region *);
  double gyration(int, double, double *);
  double gyration(int, double, double *, Region *);
  void angmom(int, double *, double *);
  void angmom(int, double *, double *, Region *);
  void torque(int, double *, double *);
  void torque(int, double *, double *, Region *);
  void inertia(int, double *, double[3][3]);
  void inertia(int, double *, double[3][3], Region *);
  void omega(double *, double[3][3], double *);
 private:
  int me;
  std::map<tagint, int> *hash;
  int find_unused();
  static void molring(int, char *, void *);
  int molbit;
};
}
#endif
