#ifdef NPAIR_CLASS
NPairStyle(half / bin / atomonly / newton, NPairHalfBinAtomonlyNewton,
           NP_HALF | NP_BIN | NP_ATOMONLY | NP_NEWTON | NP_ORTHO);
#else
#ifndef LMP_NPAIR_HALF_BIN_ATOMONLY_NEWTON_H
#define LMP_NPAIR_HALF_BIN_ATOMONLY_NEWTON_H
#include "lmptype.h"
#include "npair.h"
namespace LAMMPS_NS {
class NPairHalfBinAtomonlyNewton : public NPair {
public:
  NPairHalfBinAtomonlyNewton(class LAMMPS *);
  void build(class NeighList *) override;
};
} // namespace LAMMPS_NS
#endif
#endif
