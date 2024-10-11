#ifndef LMP_COMM_TILED_H
#define LMP_COMM_TILED_H 
#include "comm.h"
namespace LAMMPS_NS {
class CommTiled : public Comm {
 public:
  CommTiled(class LAMMPS *);
  CommTiled(class LAMMPS *, class Comm *);
  ~CommTiled() override;
  void init() override;
  void setup() override;
  void forward_comm(int dummy = 0) override;
  void reverse_comm() override;
  void exchange() override;
  void borders() override;
  void forward_comm(class Pair *) override;
  void reverse_comm(class Pair *) override;
  void forward_comm(class Fix *, int size = 0) override;
  void reverse_comm(class Fix *, int size = 0) override;
  void reverse_comm_variable(class Fix *) override;
  void forward_comm(class Compute *) override;
  void reverse_comm(class Compute *) override;
  void forward_comm_array(int, double **) override;
  int exchange_variable(int, double *, double *&) override;
  int exchange_variable_all2all(int, double *, double *&) override;
  void coord2proc_setup() override;
  int coord2proc(double *, int &, int &, int &) override;
  double memory_usage() override;
 private:
  int nswap;
  int maxswap;
  int *nsendproc, *nrecvproc;
  int *sendother, *recvother;
  int *sendself;
  int *nprocmax;
  int **sendproc, **recvproc;
  int **sendnum, **recvnum;
  int **size_forward_recv;
  int **firstrecv;
  int **size_reverse_send;
  int **size_reverse_recv;
  int **forward_recv_offset;
  int **reverse_recv_offset;
  int ***sendlist;
  int **maxsendlist;
  int **pbc_flag;
  int ***pbc;
  double ***sendbox;
  double **cutghostmulti;
  double **cutghostmultiold;
  double ****sendbox_multi;
  double ****sendbox_multiold;
  int *nexchproc;
  int *nexchprocmax;
  int **exchproc;
  int **exchnum;
  double *buf_send;
  double *buf_recv;
  int maxsend, maxrecv;
  int smaxone, rmaxone;
  int smaxall, rmaxall;
  int maxrequest;
  MPI_Request *requests;
  struct RCBinfo {
    double mysplit[3][2];
    double cutfrac;
    int dim;
  };
  RCBinfo *rcbinfo;
  int noverlap;
  int maxoverlap;
  int *overlap;
  double *prd;
  double *boxlo, *boxhi;
  double *sublo, *subhi;
  int dimension;
  void init_buffers();
  typedef void (CommTiled::*BoxDropPtr)(int, double *, double *, int &);
  BoxDropPtr box_drop;
  void box_drop_brick(int, double *, double *, int &);
  void box_drop_tiled(int, double *, double *, int &);
  void box_drop_tiled_recurse(double *, double *, int, int, int &);
  typedef void (CommTiled::*BoxOtherPtr)(int, int, int, double *, double *);
  BoxOtherPtr box_other;
  void box_other_brick(int, int, int, double *, double *);
  void box_other_tiled(int, int, int, double *, double *);
  typedef int (CommTiled::*BoxTouchPtr)(int, int, int);
  BoxTouchPtr box_touch;
  int box_touch_brick(int, int, int);
  int box_touch_tiled(int, int, int);
  typedef int (CommTiled::*PointDropPtr)(int, double *);
  PointDropPtr point_drop;
  int point_drop_brick(int, double *);
  int point_drop_tiled(int, double *);
  int point_drop_tiled_recurse(double *, int, int);
  int closer_subbox_edge(int, double *);
  void grow_send(int, int);
  void grow_recv(int);
  void grow_list(int, int, int);
  void allocate_swap(int);
  void grow_swap_send(int, int, int);
  void grow_swap_send_multi(int, int);
  void grow_swap_recv(int, int);
  void deallocate_swap(int);
};
}
#endif
