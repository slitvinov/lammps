#ifdef NPAIR_CLASS
NPairStyle(halffull/newton,
           NPairHalffullNewton,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_ORTHO | NP_TRI);
NPairStyle(halffull/newton/skip,
           NPairHalffullNewton,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_ORTHO | NP_TRI | NP_SKIP);
#else
#ifndef LMP_NPAIR_HALFFULL_NEWTON_H
#define LMP_NPAIR_HALFFULL_NEWTON_H 
#include "npair.h"
namespace LAMMPS_NS {
class NPairHalffullNewton : public NPair {
 public:
  NPairHalffullNewton(class LAMMPS *);
  void build(class NeighList *) override;
};
}
#endif
#endif
