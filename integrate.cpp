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
void Integrate::init() {
  update->atimestep = update->ntimestep;
}
void Integrate::ev_set(bigint ntimestep) {
  int i, flag;
  int tdflag = 0;
  flag = 0;
  int eflag_global = 0;
  flag = 0;
  int eflag_atom = 0;
  flag = 0;
  int vflag_global = 0;
  flag = 0;
  int vflag_atom = 0;
  flag = 0;
  int cvflag_atom = 0;
}

