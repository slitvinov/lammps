#include <map>
#include <set>
#include <unordered_set>
#include <vector>
#include <cstring>
#include <vector>
#include <cstdio>
#include <string>
#include <cmath>
#include <mpi.h>
#include "lammps.h"
#include "pointers.h"
#include "lmptype.h"
#include "verlet.h"
#include "atom.h"
#include "atom_vec.h"
#include "pointers.h"
#include "comm.h"
#include "domain.h"
#include "fix.h"
#include "fix_nve.h"
#include "force.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"
using namespace LAMMPS_NS;
Verlet::Verlet(LAMMPS *lmp, int narg, char **arg) : Pointers(lmp) {}
void Verlet::init() {
  update->atimestep = update->ntimestep;
  torqueflag = extraflag = 0;
}
void Verlet::setup(int flag) {
  update->setupflag = 1;
  atom->setup();
  domain->pbc();
  domain->reset_box();
  comm->setup();
  comm->exchange();
  if (atom->sortfreq > 0)
    atom->sort();
  comm->borders();
  neighbor->build(1);
  neighbor->ncalls = 0;
  force->setup();
  force_clear();
  force->pair->compute(0, 0);
  if (force->newton)
    comm->reverse_comm();
  update->setupflag = 0;
}
void Verlet::run(int n) {
  bigint ntimestep;
  int nflag, sortflag;
  sortflag = 1;
  for (int i = 0; i < n; i++) {
    ntimestep = ++update->ntimestep;
    lmp->fix_nve->initial_integrate(0);
    nflag = neighbor->decide();
    domain->pbc();
    comm->exchange();
    if (sortflag && ntimestep >= atom->nextsort)
      atom->sort();
    comm->borders();
    neighbor->build(1);
    force_clear();
    force->pair->compute(0, 0);
    comm->reverse_comm();
    lmp->fix_nve->final_integrate();
  }
}
void Verlet::cleanup() {
  update->update_time();
}
void Verlet::force_clear() {
  size_t nbytes;
  int nlocal = atom->nlocal;
  if (neighbor->includegroup == 0) {
    nbytes = sizeof(double) * nlocal;
    if (force->newton)
      nbytes += sizeof(double) * atom->nghost;
    if (nbytes) {
      memset(&atom->f[0][0], 0, 3 * nbytes);
    }
  }
}
