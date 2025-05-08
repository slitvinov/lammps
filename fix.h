#ifndef LMP_FIX_H
#define LMP_FIX_H
namespace LAMMPS_NS {
class Fix : protected Pointers {
  friend class Neighbor;

public:
  static int instance_total;
  char *id, *style;
  int igroup, groupbit;
  int nevery;
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
  int local_flag;
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
} // namespace LAMMPS_NS
#endif
