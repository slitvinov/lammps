#ifdef NPAIR_CLASS
NPairStyle(half/size/bin/newtoff,
           NPairHalfSizeBinNewtoff,
           NP_HALF | NP_SIZE | NP_BIN | NP_NEWTOFF | NP_ORTHO | NP_TRI);
#else
#ifndef LMP_NPAIR_HALF_SIZE_BIN_NEWTOFF_H
#define LMP_NPAIR_HALF_SIZE_BIN_NEWTOFF_H 
#include "npair.h"
namespace LAMMPS_NS {
class NPairHalfSizeBinNewtoff : public NPair {
 public:
  NPairHalfSizeBinNewtoff(class LAMMPS *);
  void build(class NeighList *) override;
};
}
#endif
#endif
