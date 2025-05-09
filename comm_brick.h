#ifndef LMP_COMM_BRICK_H
#define LMP_COMM_BRICK_H
namespace LAMMPS_NS {
class CommBrick : public Comm {
public:
  CommBrick(class LAMMPS *);
  void init() override;
  void setup() override;
  void reverse_comm() override;
  void exchange() override;
  void borders() override;

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
  virtual void grow_send(int, int);
  virtual void grow_recv(int);
  virtual void grow_list(int, int);
  virtual void allocate_swap(int);
};
} // namespace LAMMPS_NS
#endif
