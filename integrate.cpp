#include <vector>
#include <map>
#include <cstdio>
#include <mpi.h>
#include <string>
#include "lammps.h"
#include "pointers.h"
#include "lmptype.h"
#include "integrate.h"
#include "force.h"
#include "pair.h"
#include "update.h"
using namespace LAMMPS_NS;
Integrate::Integrate(LAMMPS *lmp, int, char **) : Pointers(lmp) {
}
