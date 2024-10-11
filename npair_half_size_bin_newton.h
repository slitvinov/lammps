#ifdef NPAIR_CLASS
NPairStyle(half/size/bin/newton,
           NPairHalfSizeBinNewton,
           NP_HALF | NP_SIZE | NP_BIN | NP_NEWTON | NP_ORTHO);
#else
#ifndef LMP_NPAIR_HALF_SIZE_BIN_NEWTON_H
#define LMP_NPAIR_HALF_SIZE_BIN_NEWTON_H 
#include "npair.h"
namespace LAMMPS_NS {
class NPairHalfSizeBinNewton : public NPair {
 public:
  NPairHalfSizeBinNewton(class LAMMPS *);
  void build(class NeighList *) override;
};
}
#endif
#endif
