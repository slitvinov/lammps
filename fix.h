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
  bigint next_reneighbor;
  int nevery;
  int thermo_energy;
  int thermo_virial;
  int energy_global_flag;
  int energy_peratom_flag;
  int virial_global_flag;
  int virial_peratom_flag;
  int ecouple_flag;
  int time_integrate;
  int rigid_flag;
  int no_change_box;
  int time_depend;
  int create_attribute;
  int restart_pbc;
  int wd_header;
  int wd_section;
  int dynamic_group_allow;
  int dof_flag;
  int special_alter_flag;
  int enforce2d_flag;
  int respa_level_support;
  int respa_level;
  int maxexchange;
  int maxexchange_dynamic;
  int pre_exchange_migrate;
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
  int pergrid_flag;
  int pergrid_freq;
  int extscalar;
  int extvector;
  int *extlist;
  int extarray;
  double *vector_atom;
  double **array_atom;
  double *vector_local;
  double **array_local;
  int comm_forward;
  int comm_reverse;
  int comm_border;
  double virial[6];
  double *eatom, **vatom;
  double **cvatom;
  int centroidstressflag;
  int restart_reset;
  int forward_comm_device;
  ExecutionSpace execution_space;
  unsigned int datamask_read, datamask_modify;
  Fix(class LAMMPS *, int, char **);
  virtual int setmask() = 0;
  virtual void post_constructor() {}
  virtual void init() {}
  virtual void init_list(int, class NeighList *) {}
  virtual void setup(int) {}
  virtual void setup_pre_exchange() {}
  virtual void setup_pre_neighbor() {}
  virtual void setup_post_neighbor() {}
  virtual void setup_pre_force(int) {}
  virtual void setup_pre_reverse(int, int) {}
  virtual void min_setup(int) {}
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
  virtual void restart(char *) {}
  virtual void grow_arrays(int) {}
  virtual void copy_arrays(int, int, int) {}
  virtual void set_arrays(int) {}
  virtual void update_arrays(int, int) {}
  virtual void clear_bonus() {}
  virtual int pack_border(int, int *, double *) { return 0; }
  virtual int unpack_border(int, int, double *) { return 0; }
  virtual int pack_exchange(int, double *) { return 0; }
  virtual int unpack_exchange(int, double *) { return 0; }
  virtual int pack_restart(int, double *) { return 0; }
  virtual void unpack_restart(int, int) {}
  virtual int size_restart(int) { return 0; }
  virtual int maxsize_restart() { return 0; }
  virtual void setup_pre_force_respa(int, int) {}
  virtual void initial_integrate_respa(int, int, int) {}
  virtual void post_integrate_respa(int, int) {}
  virtual void pre_force_respa(int, int, int) {}
  virtual void post_force_respa(int, int, int) {}
  virtual void final_integrate_respa(int, int) {}
  virtual void min_pre_exchange() {}
  virtual void min_pre_neighbor() {}
  virtual void min_post_neighbor() {}
  virtual void min_pre_force(int) {}
  virtual void min_pre_reverse(int, int) {}
  virtual void min_post_force(int) {}
  virtual double min_energy(double *) { return 0.0; }
  virtual void min_store() {}
  virtual void min_clearstore() {}
  virtual void min_pushstore() {}
  virtual void min_popstore() {}
  virtual int min_reset_ref() { return 0; }
  virtual void min_step(double, double *) {}
  virtual double max_alpha(double *) { return 0.0; }
  virtual int min_dof() { return 0; }
  virtual int pack_forward_comm(int, int *, double *, int, int *) { return 0; }
  virtual void unpack_forward_comm(int, int, double *) {}
  virtual int pack_reverse_comm_size(int, int) { return 0; }
  virtual int pack_reverse_comm(int, int, double *) { return 0; }
  virtual void unpack_reverse_comm(int, int *, double *) {}
  virtual void reset_grid(){};
  virtual void pack_forward_grid(int, void *, int, int *){};
  virtual void unpack_forward_grid(int, void *, int, int *){};
  virtual void pack_reverse_grid(int, void *, int, int *){};
  virtual void unpack_reverse_grid(int, void *, int, int *){};
  virtual void pack_remap_grid(int, void *, int, int *){};
  virtual void unpack_remap_grid(int, void *, int, int *){};
  virtual int unpack_read_grid(int, char *) { return 0; };
  virtual void pack_write_grid(int, void *){};
  virtual void unpack_write_grid(int, void *, int *){};
  virtual int get_grid_by_name(const std::string &, int &) { return -1; };
  virtual void *get_grid_by_index(int) { return nullptr; };
  virtual int get_griddata_by_name(int, const std::string &, int &) {
    return -1;
  };

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
