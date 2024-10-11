#ifdef NPAIR_CLASS
NPairStyle(half/bin/newtoff,
           NPairHalfBinNewtoff,
           NP_HALF | NP_BIN | NP_NEWTOFF | NP_ORTHO | NP_TRI);
#else
#ifndef LMP_NPAIR_HALF_BIN_NEWTOFF_H
#define LMP_NPAIR_HALF_BIN_NEWTOFF_H 
#include "npair.h"
namespace LAMMPS_NS {
class NPairHalfBinNewtoff : public NPair {
 public:
  NPairHalfBinNewtoff(class LAMMPS *);
  void build(class NeighList *) override;
};
}
#endif
#endif
