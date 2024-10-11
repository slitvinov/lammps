#ifndef LAMMPS_RCB_H
#define LAMMPS_RCB_H 
#include "pointers.h"
namespace LAMMPS_NS {
class RCB : protected Pointers {
 public:
  int noriginal;
  int nfinal;
  int nkeep;
  int *recvproc;
  int *recvindex;
  double *lo, *hi;
  double cut;
  int cutdim;
  int *sendproc;
  int *sendindex;
  RCB(class LAMMPS *);
  ~RCB() override;
  void compute(int, int, double **, double *, double *, double *);
  void compute_old(int, int, double **, double *, double *, double *);
  void invert(int sortflag = 0);
  double memory_usage();
  struct Median {
    double totallo, totalhi;
    double valuelo, valuehi;
    double wtlo, wthi;
    int countlo, counthi;
    int proclo, prochi;
  };
  struct BBox {
    double lo[3], hi[3];
  };
 private:
  int me, nprocs;
  struct Dot {
    double x[3];
    double wt;
    int proc;
    int index;
  };
  struct Tree {
    double cut;
    int dim;
  };
  struct Invert {
    int rindex;
    int sproc;
    int sindex;
  };
  Dot *dots;
  int ndot;
  int maxdot;
  int ndotorig;
  int nlist;
  int maxlist;
  int *dotlist;
  int *dotmark;
  int *dotmark_select;
  int maxbuf;
  Dot *buf;
  int maxrecv, maxsend;
  BBox bbox;
  MPI_Op box_op, med_op;
  MPI_Datatype box_type, med_type;
  int reuse;
  int dottop;
  double bboxlo[3];
  double bboxhi[3];
  Tree *tree;
  int counters[7];
};
}
#endif
