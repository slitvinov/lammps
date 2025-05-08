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
#include "modify.h"
#include "pair.h"
#include "update.h"
using namespace LAMMPS_NS;
Integrate::Integrate(LAMMPS *lmp, int, char **) : Pointers(lmp) {
  elist_global = elist_atom = nullptr;
  vlist_global = vlist_atom = cvlist_atom = nullptr;
}
void Integrate::init() {
  update->atimestep = update->ntimestep;
  pair_compute_flag = 1;
}
void Integrate::ev_setup() {
  delete[] elist_global;
  delete[] elist_atom;
  delete[] vlist_global;
  delete[] vlist_atom;
  delete[] cvlist_atom;
  elist_global = elist_atom = nullptr;
  vlist_global = vlist_atom = cvlist_atom = nullptr;
  nelist_global = nelist_atom = 0;
  nvlist_global = nvlist_atom = ncvlist_atom = 0;
  nelist_global = nelist_atom = 0;
  nvlist_global = nvlist_atom = ncvlist_atom = 0;
}
void Integrate::ev_set(bigint ntimestep) {
  int i, flag;
  int tdflag = 0;
  flag = 0;
  int eflag_global = 0;
  flag = 0;
  int eflag_atom = 0;
  eflag = eflag_global + eflag_atom;
  flag = 0;
  int vflag_global = 0;
  flag = 0;
  int vflag_atom = 0;
  flag = 0;
  int cvflag_atom = 0;
  vflag = vflag_global + vflag_atom + cvflag_atom;
}

