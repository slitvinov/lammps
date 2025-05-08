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
#include "integrate.h"
#include "verlet.h"
#include "atom.h"
#include "atom_vec.h"
#include "pointers.h"
#include "comm.h"
#include "domain.h"
#include "fix.h"
#include "force.h"
#include "modify.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"
using namespace LAMMPS_NS;
Verlet::Verlet(LAMMPS *lmp, int narg, char **arg) : Integrate(lmp, narg, arg) {}
void Verlet::init() {
  Integrate::init();
  bool do_time_integrate = false;
  for (const auto &fix : modify->get_fix_list())
    if (fix->time_integrate)
      do_time_integrate = true;
  torqueflag = extraflag = 0;
}
void Verlet::setup(int flag) {
  update->setupflag = 1;
  atom->setup();
  modify->setup_pre_exchange();
  domain->pbc();
  domain->reset_box();
  comm->setup();
  comm->exchange();
  if (atom->sortfreq > 0)
    atom->sort();
  comm->borders();
  modify->setup_pre_neighbor();
  neighbor->build(1);
  modify->setup_post_neighbor();
  neighbor->ncalls = 0;
  force->setup();
  ev_set(update->ntimestep);
  force_clear();
  force->pair->compute(0, 0);
  if (force->newton)
    comm->reverse_comm();
  modify->setup(0);
  update->setupflag = 0;
}
void Verlet::run(int n) {
  bigint ntimestep;
  int nflag, sortflag;
  int n_post_integrate = modify->n_post_integrate;
  int n_pre_exchange = modify->n_pre_exchange;
  int n_pre_neighbor = modify->n_pre_neighbor;
  int n_post_neighbor = modify->n_post_neighbor;
  int n_pre_force = modify->n_pre_force;
  int n_pre_reverse = modify->n_pre_reverse;
  int n_post_force_any = modify->n_post_force_any;
  int n_end_of_step = modify->n_end_of_step;
  sortflag = 1;
  for (int i = 0; i < n; i++) {
    ntimestep = ++update->ntimestep;
    ev_set(ntimestep);
    modify->initial_integrate(0);
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
    modify->final_integrate();
  }
}
void Verlet::cleanup() {
  modify->post_run();
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
