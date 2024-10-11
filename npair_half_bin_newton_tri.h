#ifdef NPAIR_CLASS
NPairStyle(half/bin/newton/tri,
           NPairHalfBinNewtonTri,
           NP_HALF | NP_BIN | NP_NEWTON | NP_TRI);
#else
#ifndef LMP_NPAIR_HALF_BIN_NEWTON_TRI_H
#define LMP_NPAIR_HALF_BIN_NEWTON_TRI_H 
#include "npair.h"
namespace LAMMPS_NS {
class NPairHalfBinNewtonTri : public NPair {
 public:
  NPairHalfBinNewtonTri(class LAMMPS *);
  void build(class NeighList *) override;
};
}
#endif
#endif
