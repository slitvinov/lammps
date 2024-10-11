#ifndef LMP_COMM_H
#define LMP_COMM_H 
#include "pointers.h"
namespace LAMMPS_NS {
class Comm : protected Pointers {
 public:
  enum { BRICK, TILED };
  int style;
  enum { LAYOUT_UNIFORM, LAYOUT_NONUNIFORM, LAYOUT_TILED };
  int layout;
  enum { SINGLE, MULTI, MULTIOLD };
  int mode;
  int me, nprocs;
  int ghost_velocity;
  double cutghost[3];
  double cutghostuser;
  double *cutusermulti;
  double *cutusermultiold;
  int ncollections;
  int ncollections_cutoff;
  int recv_from_partition;
  int send_to_partition;
  int other_partition_style;
  int nthreads;
  int procgrid[3];
  int user_procgrid[3];
  int myloc[3];
  int procneigh[3][2];
  double *xsplit, *ysplit, *zsplit;
  int ***grid2proc;
  int rcbnew;
  double mysplit[3][2];
  double rcbcutfrac;
  int rcbcutdim;
  Comm(class LAMMPS *);
  ~Comm() override;
  void copy_arrays(class Comm *);
  virtual void init();
  void modify_params(int, char **);
  void set_processors(int, char **);
  virtual void set_proc_grid(int outflag = 1);
  double get_comm_cutoff();
  virtual void setup() = 0;
  virtual void forward_comm(int dummy = 0) = 0;
  virtual void reverse_comm() = 0;
  virtual void exchange() = 0;
  virtual void borders() = 0;
  virtual void forward_comm(class Pair *) = 0;
  virtual void reverse_comm(class Pair *) = 0;
  virtual void forward_comm(class Fix *, int size = 0) = 0;
  virtual void reverse_comm(class Fix *, int size = 0) = 0;
  virtual void reverse_comm_variable(class Fix *) = 0;
  virtual void forward_comm(class Compute *) = 0;
  virtual void reverse_comm(class Compute *) = 0;
  virtual void forward_comm_array(int, double **) = 0;
  virtual int exchange_variable(int, double *, double *&) = 0;
  virtual int exchange_variable_all2all(int, double *, double *&) = 0;
  virtual void coord2proc_setup() {}
  virtual int coord2proc(double *, int &, int &, int &);
  virtual double memory_usage() = 0;
  void ring(int, int, void *, int, void (*)(int, char *, void *), void *, void *, int self = 1);
  int rendezvous(int, int, char *, int, int, int *,
                 int (*)(int, char *, int &, int *&, char *&, void *), int, char *&, int, void *,
                 int statflag = 0);
  virtual void *extract(const char *, int &) { return nullptr; }
 protected:
  int bordergroup;
  int triclinic;
  int map_style;
  int comm_x_only, comm_f_only;
  int size_forward;
  int size_reverse;
  int size_border;
  int maxforward, maxreverse;
  int maxexchange;
  int maxexchange_atom;
  int maxexchange_fix;
  int maxexchange_fix_dynamic;
  int bufextra;
  int gridflag;
  int mapflag;
  char xyz[4];
  char *customfile;
  char *outfile;
  int otherflag;
  int other_style;
  int other_procgrid[3];
  int other_coregrid[3];
  int ncores;
  int coregrid[3];
  int user_coregrid[3];
  int multi_reduce;
  void init_exchange();
  int rendezvous_all2all(int, char *, int, int, int *,
                         int (*)(int, char *, int &, int *&, char *&, void *), int, char *&, int,
                         void *, int);
  void rendezvous_stats(int, int, int, int, int, int, bigint);
 public:
  enum { MULTIPLE };
};
}
#endif
