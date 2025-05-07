#ifndef LMP_ATOM_VEC_H
#define LMP_ATOM_VEC_H
namespace LAMMPS_NS {
class AtomVec : protected Pointers {
public:
  enum { PER_ATOM = 0, PER_TYPE = 1 };
  int molecular;
  int bonds_allow;
  int mass_type;
  int dipole_type;
  int forceclearflag;
  int comm_x_only;
  int comm_f_only;
  int size_forward;
  int size_reverse;
  int size_border;
  int size_velocity;
  int size_data_atom;
  int size_data_vel;
  int xcol_data;
  int maxexchange;
  int bonus_flag;
  int size_forward_bonus;
  int size_border_bonus;
  int size_restart_bonus_one;
  int size_data_bonus;
  int nargcopy;
  char **argcopy;
  std::vector<std::string> fields_grow, fields_copy, fields_comm,
      fields_comm_vel;
  std::vector<std::string> fields_reverse, fields_border, fields_border_vel;
  std::vector<std::string> fields_exchange, fields_restart, fields_create;
  std::vector<std::string> fields_data_atom, fields_data_vel;
  AtomVec(class LAMMPS *);
  void store_args(int, char **);
  virtual void process_args(int, char **);
  virtual void init();
  virtual void grow(int);
  virtual void grow_pointers() {}
  virtual void copy(int, int, int);
  virtual void clear_bonus() {}
  virtual void unpack_reverse(int, int *, double *);
  virtual int pack_border_vel(int, int *, double *, int, int *);
  virtual void unpack_border_vel(int, int, double *);
  virtual void create_atom(int, double *);
  virtual void create_atom_post(int) {}

protected:
  int nmax;
  int deform_vremap;
  int deform_groupbit;
  double *h_rate;
  tagint *tag;
  int *type, *mask;
  imageint *image;
  double **x, **v, **f;
  static const std::vector<std::string> default_grow, default_copy,
      default_comm, default_comm_vel;
  static const std::vector<std::string> default_reverse, default_border,
      default_border_vel;
  static const std::vector<std::string> default_exchange, default_restart,
      default_create;
  static const std::vector<std::string> default_data_atom, default_data_vel;
  struct Method {
    std::vector<void *> pdata;
    std::vector<int> datatype;
    std::vector<int> cols;
    std::vector<int *> maxcols;
    std::vector<int> collength;
    std::vector<void *> plength;
    std::vector<int> index;
    void resize(int nfield);
  };
  Method mgrow, mcopy;
  Method mcomm, mcomm_vel, mreverse, mborder, mborder_vel, mexchange, mrestart;
  Method mcreate, mdata_atom, mdata_vel;
  int ngrow, ncopy;
  int ncomm, ncomm_vel, nreverse, nborder, nborder_vel, nexchange, nrestart;
  int ncreate, ndata_atom, ndata_vel;
  bool *threads;
  void grow_nmax();
  void setup_fields();
  int process_fields(const std::vector<std::string> &,
                     const std::vector<std::string> &, Method *);
  void init_method(int, Method *);
};
} // namespace LAMMPS_NS
#endif
