#ifndef LMP_COMM_BRICK_H
#define LMP_COMM_BRICK_H 
#include "comm.h"
namespace LAMMPS_NS {
class CommBrick : public Comm {
 public:
  CommBrick(class LAMMPS *);
  CommBrick(class LAMMPS *, class Comm *);
  ~CommBrick() override;
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
  void *extract(const char *, int &) override;
  double memory_usage() override;
 protected:
  int nswap;
  int recvneed[3][2];
  int sendneed[3][2];
  int maxneed[3];
  int maxswap;
  int *sendnum, *recvnum;
  int *sendproc, *recvproc;
  int *size_forward_recv;
  int *size_reverse_send;
  int *size_reverse_recv;
  double *slablo, *slabhi;
  double **multilo, **multihi;
  double **multioldlo, **multioldhi;
  double **cutghostmulti;
  double **cutghostmultiold;
  int *pbc_flag;
  int **pbc;
  int *firstrecv;
  int **sendlist;
  int *localsendlist;
  int *maxsendlist;
  double *buf_send;
  double *buf_recv;
  int maxsend, maxrecv;
  int smax, rmax;
  void init_buffers();
  int updown(int, int, int, double, int, double *);
  virtual void grow_send(int, int);
  virtual void grow_recv(int);
  virtual void grow_list(int, int);
  virtual void grow_swap(int);
  virtual void allocate_swap(int);
  virtual void allocate_multi(int);
  virtual void allocate_multiold(int);
  virtual void free_swap();
  virtual void free_multi();
  virtual void free_multiold();
};
}
#endif
