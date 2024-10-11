#ifdef FIX_CLASS
FixStyle(nph/sphere,FixNPHSphere);
#else
#ifndef LMP_FIX_NPH_SPHERE_H
#define LMP_FIX_NPH_SPHERE_H 
#include "fix_nh_sphere.h"
namespace LAMMPS_NS {
class FixNPHSphere : public FixNHSphere {
 public:
  FixNPHSphere(class LAMMPS *, int, char **);
};
}
#endif
#endif
