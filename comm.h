#ifndef LMP_COMM_H
#define LMP_COMM_H
namespace LAMMPS_NS {
class Comm : protected Pointers {
public:
  enum { BRICK };
  int style;
  enum { LAYOUT_UNIFORM };
  int layout;
  enum { SINGLE };
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
  virtual void init();
  void modify_params(int, char **);
  virtual void set_proc_grid(int outflag = 1);
  double get_comm_cutoff();
  virtual void setup() = 0;
  virtual void reverse_comm() = 0;
  virtual void exchange() = 0;
  virtual void borders() = 0;
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

public:
  enum { MULTIPLE };
};
} // namespace LAMMPS_NS
#endif
