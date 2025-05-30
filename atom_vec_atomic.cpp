#include <map>
#include <set>
#include <vector>
#include <cstdio>
#include <mpi.h>
#include <string>
#include "lammps.h"
#include "pointers.h"
#include "lmptype.h"
#include "atom_vec.h"
#include "atom_vec_atomic.h"
#include "atom.h"
#include "pointers.h"
using namespace LAMMPS_NS;
AtomVecAtomic::AtomVecAtomic(LAMMPS *lmp) : AtomVec(lmp) {
  molecular = Atom::ATOMIC;
  mass_type = PER_TYPE;
  fields_data_atom = {"id", "type", "x"};
  fields_data_vel = {"id", "v"};
  setup_fields();
}
