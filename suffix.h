#ifndef LMP_SUFFIX_H
#define LMP_SUFFIX_H 
namespace LAMMPS_NS {
namespace Suffix {
  enum { NONE = 0, OPT = 1 << 0, GPU = 1 << 1, OMP = 1 << 2, INTEL = 1 << 3, KOKKOS = 1 << 4 };
}
}
#endif
