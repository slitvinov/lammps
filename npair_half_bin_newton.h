#ifdef NPAIR_CLASS
NPairStyle(half/bin/newton,
           NPairHalfBinNewton,
           NP_HALF | NP_BIN | NP_MOLONLY | NP_NEWTON | NP_ORTHO);
#else
#ifndef LMP_NPAIR_HALF_BIN_NEWTON_H
#define LMP_NPAIR_HALF_BIN_NEWTON_H 
#include "npair.h"
namespace LAMMPS_NS {
class NPairHalfBinNewton : public NPair {
 public:
  NPairHalfBinNewton(class LAMMPS *);
  void build(class NeighList *) override;
};
}
#endif
#endif
