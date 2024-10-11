#ifdef NPAIR_CLASS
NPairStyle(full/bin/atomonly,
           NPairFullBinAtomonly,
           NP_FULL | NP_BIN | NP_ATOMONLY |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);
#else
#ifndef LMP_NPAIR_FULL_BIN_ATOMONLY_H
#define LMP_NPAIR_FULL_BIN_ATOMONLY_H 
#include "npair.h"
namespace LAMMPS_NS {
class NPairFullBinAtomonly : public NPair {
 public:
  NPairFullBinAtomonly(class LAMMPS *);
  void build(class NeighList *) override;
};
}
#endif
#endif
