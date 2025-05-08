#ifndef LMP_FIX_H
#define LMP_FIX_H
namespace LAMMPS_NS {
class Fix : protected Pointers {
  friend class Neighbor;

public:
  static int instance_total;
  char *id, *style;
  int igroup, groupbit;
  int restart_global;
  int restart_peratom;
  int restart_file;
  int force_reneighbor;
  int box_change;
  enum {
    NO_BOX_CHANGE = 0,
    BOX_CHANGE_ANY = 1 << 0,
    BOX_CHANGE_DOMAIN = 1 << 1,
    BOX_CHANGE_X = 1 << 2,
    BOX_CHANGE_Y = 1 << 3,
    BOX_CHANGE_Z = 1 << 4,
    BOX_CHANGE_YZ = 1 << 5,
    BOX_CHANGE_XZ = 1 << 6,
    BOX_CHANGE_XY = 1 << 7,
    BOX_CHANGE_SIZE = BOX_CHANGE_X | BOX_CHANGE_Y | BOX_CHANGE_Z,
    BOX_CHANGE_SHAPE = BOX_CHANGE_YZ | BOX_CHANGE_XZ | BOX_CHANGE_XY
  };
  int nevery;
  int time_integrate;
  int maxexchange;
  int stores_ids;
  int scalar_flag;
  int vector_flag;
  int array_flag;
  int size_vector;
  int size_array_rows;
  int size_array_cols;
  int size_vector_variable;
  int size_array_rows_variable;
  int global_freq;
  int peratom_flag;
  int size_peratom_cols;
  int peratom_freq;
  int local_flag;
  int size_local_rows;
  int size_local_cols;
  int local_freq;
  int extscalar;
  int extvector;
  int *extlist;
  int extarray;
  double *vector_atom;
  double **array_atom;
  int comm_forward;
  int comm_reverse;
  int comm_border;
  double virial[6];
  double *eatom, **vatom;
  double **cvatom;
  unsigned int datamask_read, datamask_modify;
  Fix(class LAMMPS *, int, char **);
  virtual void post_constructor() {}
  virtual void init() {}
  virtual void init_list(int, class NeighList *) {}
  virtual void setup(int) {}
  virtual void initial_integrate(int) {}
  virtual void post_integrate() {}
  virtual void pre_exchange() {}
  virtual void pre_neighbor() {}
  virtual void post_neighbor() {}
  virtual void pre_force(int) {}
  virtual void pre_reverse(int, int) {}
  virtual void post_force(int) {}
  virtual void final_integrate() {}
  virtual void post_run() {}

protected:
  int instance_me;
  int evflag;
  int eflag_either, eflag_global, eflag_atom;
  int vflag_either, vflag_global, vflag_atom, cvflag_atom;
  int maxeatom, maxvatom, maxcvatom;
  int copymode;
  int dynamic;
};
namespace FixConst {
enum {
  INITIAL_INTEGRATE = 1 << 0,
  POST_INTEGRATE = 1 << 1,
  PRE_EXCHANGE = 1 << 2,
  PRE_NEIGHBOR = 1 << 3,
  POST_NEIGHBOR = 1 << 4,
  PRE_FORCE = 1 << 5,
  PRE_REVERSE = 1 << 6,
  POST_FORCE = 1 << 7,
  FINAL_INTEGRATE = 1 << 8,
  END_OF_STEP = 1 << 9,
  POST_RUN = 1 << 10,
  INITIAL_INTEGRATE_RESPA = 1 << 11,
  POST_INTEGRATE_RESPA = 1 << 12,
  PRE_FORCE_RESPA = 1 << 13,
  POST_FORCE_RESPA = 1 << 14,
  FINAL_INTEGRATE_RESPA = 1 << 15,
  MIN_PRE_EXCHANGE = 1 << 16,
  MIN_PRE_NEIGHBOR = 1 << 17,
  MIN_POST_NEIGHBOR = 1 << 18,
  MIN_PRE_FORCE = 1 << 19,
  MIN_PRE_REVERSE = 1 << 20,
  MIN_POST_FORCE = 1 << 21,
  MIN_ENERGY = 1 << 22
};
}
} // namespace LAMMPS_NS
#endif
