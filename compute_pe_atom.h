#ifdef COMPUTE_CLASS
ComputeStyle(pe/atom,ComputePEAtom);
#else
#ifndef LMP_COMPUTE_PE_ATOM_H
#define LMP_COMPUTE_PE_ATOM_H 
#include "compute.h"
namespace LAMMPS_NS {
class ComputePEAtom : public Compute {
 public:
  ComputePEAtom(class LAMMPS *, int, char **);
  ~ComputePEAtom() override;
  void init() override {}
  void compute_peratom() override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double memory_usage() override;
 private:
  int pairflag;
  int fixflag;
  int nmax;
  double *energy;
};
}
#endif
#endif
