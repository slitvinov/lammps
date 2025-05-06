#ifndef LMP_COMPUTE_H
#define LMP_COMPUTE_H
namespace LAMMPS_NS {
class Compute : protected Pointers {
  friend class Neighbor;

public:
  enum {
    INVOKED_NONE = 0,
    INVOKED_SCALAR = 1 << 0,
    INVOKED_VECTOR = 1 << 1,
    INVOKED_ARRAY = 1 << 2,
    INVOKED_PERATOM = 1 << 3,
    INVOKED_LOCAL = 1 << 4,
    INVOKED_PERGRID = 1 << 5,
  };
  static int instance_total;
  char *id, *style;
  int igroup, groupbit;
  double scalar;
  double *vector;
  double **array;
  double *vector_atom;
  double **array_atom;
  double *vector_local;
  double **array_local;
  int scalar_flag;
  int vector_flag;
  int array_flag;
  int size_vector;
  int size_array_rows;
  int size_array_cols;
  int size_vector_variable;
  int size_array_rows_variable;
  int peratom_flag;
  int size_peratom_cols;
  int local_flag;
  int size_local_rows;
  int size_local_cols;
  int pergrid_flag;
  int extscalar;
  int extvector;
  int *extlist;
  int extarray;
  int tempflag;
  int pressflag;
  int pressatomflag;
  int peflag;
  int peatomflag;
  int create_attribute;
  int tempbias;
  int timeflag;
  int ntime;
  int maxtime;
  bigint *tlist;
  int invoked_flag;
  bigint invoked_scalar;
  bigint invoked_vector;
  bigint invoked_array;
  bigint invoked_peratom;
  bigint invoked_local;
  bigint invoked_pergrid;
  double dof;
  int comm_forward;
  int comm_reverse;
  int dynamic_group_allow;
  ExecutionSpace execution_space;
  unsigned int datamask_read, datamask_modify;
  int copymode;
  Compute(class LAMMPS *, int, char **);
  ~Compute() override;
  void modify_params(int, char **);
  virtual void reset_extra_dof();
  virtual void init() = 0;
  virtual void init_list(int, class NeighList *) {}
  virtual void setup() {}
  virtual double compute_scalar() { return 0.0; }
  virtual void compute_vector() {}
  virtual void compute_array() {}
  virtual void compute_peratom() {}
  virtual void compute_local() {}
  virtual void compute_pergrid() {}
  virtual void set_arrays(int) {}
  virtual int pack_forward_comm(int, int *, double *, int, int *) { return 0; }
  virtual void unpack_forward_comm(int, int, double *) {}
  virtual int pack_reverse_comm(int, int, double *) { return 0; }
  virtual void unpack_reverse_comm(int, int *, double *) {}
  virtual void reset_grid(){};
  virtual int get_grid_by_name(const std::string &, int &) { return -1; };
  virtual void *get_grid_by_index(int) { return nullptr; };
  virtual int get_griddata_by_name(int, const std::string &, int &) {
    return -1;
  };
  virtual void *get_griddata_by_index(int) { return nullptr; };
  virtual void dof_remove_pre() {}
  virtual int dof_remove(int) { return 0; }
  virtual void remove_bias(int, double *) {}
  virtual void remove_bias_thr(int, double *, double *) {}
  virtual void remove_bias_all() {}
  virtual void reapply_bias_all() {}
  virtual void restore_bias(int, double *) {}
  virtual void restore_bias_thr(int, double *, double *) {}
  virtual void restore_bias_all() {}
  virtual void reset_extra_compute_fix(const char *);
  virtual void lock_enable() {}
  virtual void lock_disable() {}
  virtual int lock_length() { return 0; }
  virtual void lock(class Fix *, bigint, bigint) {}
  virtual void unlock(class Fix *) {}
  virtual void refresh() {}
  void addstep(bigint);
  int matchstep(bigint);
  void clearstep();
  virtual void pair_setup_callback(int, int) {}
  virtual void pair_tally_callback(int, int, int, int, double, double, double,
                                   double, double, double) {}

protected:
  int instance_me;
  double natoms_temp;
  double extra_dof;
  int fix_dof;
  int dynamic;
  int dynamic_user;
  double vbias[3];
  double **vbiasall;
  int maxbias;
  inline int sbmask(int j) const { return j >> SBBITS & 3; }
  void adjust_dof_fix();
};
} // namespace LAMMPS_NS
#endif
