#ifdef NPAIR_CLASS
NPairStyle(halffull/newtoff,
           NPairHalffullNewtoff,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI);
NPairStyle(halffull/newtoff/skip,
           NPairHalffullNewtoff,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_SKIP);
NPairStyle(halffull/newtoff/ghost,
           NPairHalffullNewtoff,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_GHOST);
NPairStyle(halffull/newtoff/skip/ghost,
           NPairHalffullNewtoff,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_GHOST);
#else
#ifndef LMP_NPAIR_HALFFULL_NEWTOFF_H
#define LMP_NPAIR_HALFFULL_NEWTOFF_H 
#include "npair.h"
namespace LAMMPS_NS {
class NPairHalffullNewtoff : public NPair {
 public:
  NPairHalffullNewtoff(class LAMMPS *);
  void build(class NeighList *) override;
};
}
#endif
#endif
