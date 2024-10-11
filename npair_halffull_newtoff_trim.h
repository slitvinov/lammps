#ifdef NPAIR_CLASS
NPairStyle(halffull/newtoff/trim,
           NPairHalffullNewtoffTrim,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_TRIM);
NPairStyle(halffull/newtoff/skip/trim,
           NPairHalffullNewtoffTrim,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_TRIM);
NPairStyle(halffull/newtoff/ghost/trim,
           NPairHalffullNewtoffTrim,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_TRIM);
NPairStyle(halffull/newtoff/skip/ghost/trim,
           NPairHalffullNewtoffTrim,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_GHOST | NP_TRIM);
#else
#ifndef LMP_NPAIR_HALFFULL_NEWTOFF_TRIM_H
#define LMP_NPAIR_HALFFULL_NEWTOFF_TRIM_H 
#include "npair.h"
namespace LAMMPS_NS {
class NPairHalffullNewtoffTrim : public NPair {
 public:
  NPairHalffullNewtoffTrim(class LAMMPS *);
  void build(class NeighList *) override;
};
}
#endif
#endif
