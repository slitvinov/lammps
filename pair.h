#ifndef LMP_PAIR_H
#define LMP_PAIR_H
namespace LAMMPS_NS {
class Pair : protected Pointers {
  friend class FixQEq;
  friend class PairHybrid;
  friend class PairHybridScaled;
  friend class Info;
  friend class Neighbor;

public:
  static int instance_total;
  double eng_vdwl, eng_coul;
  double virial[6];
  double *eatom, **vatom;
  double **cvatom;
  double cutforce;
  double **cutsq;
  int **setflag;
  int comm_forward;
  int comm_reverse;
  int comm_reverse_off;
  int single_enable;
  int born_matrix_enable;
  int restartinfo;
  int one_coeff;
  int manybody_flag;
  int unit_convert_flag;
  int no_virial_fdotr_compute;
  int writedata;
  int finitecutflag;
  int ghostneigh;
  double **cutghost;
  int ewaldflag;
  int msmflag;
  int tip4pflag;
  int spinflag;
  int reinitflag;
  int tail_flag;
  double etail, ptail;
  double etail_ij, ptail_ij;
  int trim_flag;
  int evflag;
  int eflag_either, eflag_global, eflag_atom;
  int vflag_either, vflag_global, vflag_atom, cvflag_atom;
  int ncoultablebits;
  int ndisptablebits;
  double tabinnersq;
  double tabinnerdispsq;
  double *rtable, *drtable, *ftable, *dftable, *ctable, *dctable;
  double *etable, *detable, *ptable, *dptable, *vtable, *dvtable;
  double *rdisptable, *drdisptable, *fdisptable, *dfdisptable;
  double *edisptable, *dedisptable;
  int ncoulshiftbits, ncoulmask;
  int ndispshiftbits, ndispmask;
  int nextra;
  double *pvector;
  int single_extra;
  double *svector;
  class NeighList *list;
  class NeighList *listhalf;
  class NeighList *listfull;
  int allocated;
  int compute_flag;
  int mixed_flag;
  bool did_mix;
  enum { GEOMETRIC, ARITHMETIC, SIXTHPOWER };
  int beyond_contact, nondefault_history_transfer;
  unsigned int datamask_read, datamask_modify;
  Pair(class LAMMPS *);
  void init();
  virtual void setup() {}
  void ev_tally(int, int, int, int, double, double, double, double, double,
                double);
  void ev_tally3(int, int, int, double, double, double *, double *, double *,
                 double *);
  void v_tally2_newton(int, double *, double *);
  void v_tally3(int, int, int, double *, double *, double *, double *);
  void v_tally4(int, int, int, int, double *, double *, double *, double *,
                double *, double *);
  virtual void compute(int, int) = 0;
  virtual void compute_inner() {}
  virtual void compute_middle() {}
  virtual void compute_outer(int, int) {}
  virtual void born_matrix(int, int, int, int, double, double, double,
                           double &du, double &du2) {
    du = du2 = 0.0;
  }
  virtual void finish() {}
  virtual void settings(int, char **) = 0;
  virtual void coeff(int, char **) = 0;
  virtual void init_style();
  virtual void init_list(int, class NeighList *);
  virtual double init_one(int, int) { return 0.0; }
  virtual int pack_forward_comm(int, int *, double *, int, int *) { return 0; }
  virtual void unpack_forward_comm(int, int, double *) {}
  virtual int pack_reverse_comm(int, int, double *) { return 0; }
  virtual void unpack_reverse_comm(int, int *, double *) {}
  virtual void reset_grid() {}
  virtual void pack_forward_grid(int, void *, int, int *) {}
  virtual void unpack_forward_grid(int, void *, int, int *) {}
  virtual void pack_reverse_grid(int, void *, int, int *) {}
  virtual void unpack_reverse_grid(int, void *, int, int *) {}
  void set_copymode(int value) { copymode = value; }
  virtual void *extract(const char *, int &) { return nullptr; }
  virtual void *extract_peratom(const char *, int &) { return nullptr; }
  virtual void swap_eam(double *, double **) {}
  virtual void min_xf_pointers(int, double **, double **) {}
  virtual void min_xf_get(int) {}
  virtual void min_x_set(int) {}
  virtual void transfer_history(double *, double *, int, int) {}
  virtual double atom2cut(int) { return 0.0; }
  virtual double radii2cut(double, double) { return 0.0; }

protected:
  int num_tally_compute;
  class Compute **list_tally_compute;

public:
protected:
  int instance_me;
  int offset_flag, mix_flag;
  double tabinner;
  double tabinner_disp;

protected:
  int nelements;
  char **elements;
  int *elem1param;
  int **elem2param;
  int ***elem3param;
  int *map;
  int nparams;
  int maxparam;

public:
  typedef union {
    int i;
    float f;
  } union_int_float_t;
  inline int fdotr_is_set() const { return vflag_fdotr; }

protected:
  int vflag_fdotr;
  int maxeatom, maxvatom, maxcvatom;
  int copymode;
  inline int sbmask(int j) const { return j >> SBBITS & 3; }
};
} // namespace LAMMPS_NS
#endif
