#ifdef NPAIR_CLASS
NPairStyle(copy,
           NPairCopy,
           NP_COPY);
#else
#ifndef LMP_NPAIR_COPY_H
#define LMP_NPAIR_COPY_H 
#include "npair.h"
namespace LAMMPS_NS {
class NPairCopy : public NPair {
 public:
  NPairCopy(class LAMMPS *);
  void build(class NeighList *) override;
};
}
#endif
#endif
