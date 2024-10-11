#ifdef NPAIR_CLASS
NPairStyle(full/bin,
           NPairFullBin,
           NP_FULL | NP_BIN | NP_MOLONLY |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);
#else
#ifndef LMP_NPAIR_FULL_BIN_H
#define LMP_NPAIR_FULL_BIN_H 
#include "npair.h"
namespace LAMMPS_NS {
class NPairFullBin : public NPair {
 public:
  NPairFullBin(class LAMMPS *);
  void build(class NeighList *) override;
};
}
#endif
#endif
