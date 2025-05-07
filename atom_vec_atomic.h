#ifndef LMP_ATOM_VEC_ATOMIC_H
#define LMP_ATOM_VEC_ATOMIC_H
namespace LAMMPS_NS {
class AtomVecAtomic : virtual public AtomVec {
public:
  AtomVecAtomic(class LAMMPS *);
};
} // namespace LAMMPS_NS
#endif
