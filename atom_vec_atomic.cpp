#include "atom_vec_atomic.h"
#include "atom.h"
using namespace LAMMPS_NS;
AtomVecAtomic::AtomVecAtomic(LAMMPS *lmp) : AtomVec(lmp) {
  molecular = Atom::ATOMIC;
  mass_type = PER_TYPE;
  fields_data_atom = {"id", "type", "x"};
  fields_data_vel = {"id", "v"};
  setup_fields();
}
