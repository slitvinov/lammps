#ifdef NPAIR_CLASS
NPairStyle(halffull/newton/trim,
           NPairHalffullNewtonTrim,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_ORTHO | NP_TRI | NP_TRIM);
NPairStyle(halffull/newton/skip/trim,
           NPairHalffullNewtonTrim,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_TRIM);
#else
#ifndef LMP_NPAIR_HALFFULL_NEWTON_TRIM_H
#define LMP_NPAIR_HALFFULL_NEWTON_TRIM_H 
#include "npair.h"
namespace LAMMPS_NS {
class NPairHalffullNewtonTrim : public NPair {
 public:
  NPairHalffullNewtonTrim(class LAMMPS *);
  void build(class NeighList *) override;
};
}
#endif
#endif
