#ifdef NPAIR_CLASS
NPairStyle(half/size/bin/newton/tri,
           NPairHalfSizeBinNewtonTri,
           NP_HALF | NP_SIZE | NP_BIN | NP_NEWTON | NP_TRI);
#else
#ifndef LMP_NPAIR_HALF_SIZE_BIN_NEWTON_TRI_H
#define LMP_NPAIR_HALF_SIZE_BIN_NEWTON_TRI_H 
#include "npair.h"
namespace LAMMPS_NS {
class NPairHalfSizeBinNewtonTri : public NPair {
 public:
  NPairHalfSizeBinNewtonTri(class LAMMPS *);
  void build(class NeighList *) override;
};
}
#endif
#endif
