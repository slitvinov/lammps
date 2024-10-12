#include "verlet.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "modify.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"
#include <cstring>
using namespace LAMMPS_NS;
Verlet::Verlet(LAMMPS *lmp, int narg, char **arg) :
  Integrate(lmp, narg, arg) {}
void Verlet::init()
{
  Integrate::init();
  bool do_time_integrate = false;
  for (const auto &fix : modify->get_fix_list())
    if (fix->time_integrate) do_time_integrate = true;
  if (!do_time_integrate && (comm->me == 0))
    error->warning(FLERR,"No fixes with time integration, atoms won't move");
  if (force->newton_pair) virial_style = VIRIAL_FDOTR;
  else virial_style = VIRIAL_PAIR;
  ev_setup();
  if (modify->get_fix_by_id("package_omp")) external_force_clear = 1;
  torqueflag = extraflag = 0;
  if (atom->torque_flag) torqueflag = 1;
  if (atom->avec->forceclearflag) extraflag = 1;
  triclinic = domain->triclinic;
}
void Verlet::setup(int flag)
{
  if (comm->me == 0 && screen) {
    fputs("Setting up Verlet run ...\n",screen);
  }
  update->setupflag = 1;
  atom->setup();
  modify->setup_pre_exchange();
  if (triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  comm->setup();
  if (neighbor->style) neighbor->setup_bins();
  comm->exchange();
  if (atom->sortfreq > 0) atom->sort();
  comm->borders();
  if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
  modify->setup_pre_neighbor();
  neighbor->build(1);
  modify->setup_post_neighbor();
  neighbor->ncalls = 0;
  force->setup();
  ev_set(update->ntimestep);
  force_clear();
  modify->setup_pre_force(vflag);
  if (pair_compute_flag) force->pair->compute(eflag,vflag);
  else if (force->pair) force->pair->compute_dummy(eflag,vflag);
  modify->setup_pre_reverse(eflag,vflag);
  if (force->newton) comm->reverse_comm();
  modify->setup(vflag);
  update->setupflag = 0;
}
void Verlet::run(int n)
{
  bigint ntimestep;
  int nflag,sortflag;
  int n_post_integrate = modify->n_post_integrate;
  int n_pre_exchange = modify->n_pre_exchange;
  int n_pre_neighbor = modify->n_pre_neighbor;
  int n_post_neighbor = modify->n_post_neighbor;
  int n_pre_force = modify->n_pre_force;
  int n_pre_reverse = modify->n_pre_reverse;
  int n_post_force_any = modify->n_post_force_any;
  int n_end_of_step = modify->n_end_of_step;
  if (atom->sortfreq > 0) sortflag = 1;
  else sortflag = 0;
  for (int i = 0; i < n; i++) {
    ntimestep = ++update->ntimestep;
    ev_set(ntimestep);
    modify->initial_integrate(vflag);
    if (n_post_integrate) modify->post_integrate();
    nflag = neighbor->decide();
    if (nflag == 0) {
      comm->forward_comm();
    } else {
      if (n_pre_exchange) {
        modify->pre_exchange();
      }
      if (triclinic) domain->x2lamda(atom->nlocal);
      domain->pbc();
      comm->exchange();
      if (sortflag && ntimestep >= atom->nextsort) atom->sort();
      comm->borders();
      if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
      if (n_pre_neighbor) {
        modify->pre_neighbor();
      }
      neighbor->build(1);
      if (n_post_neighbor) {
        modify->post_neighbor();
      }
    }
    force_clear();
    if (n_pre_force) {
      modify->pre_force(vflag);
    }
    if (pair_compute_flag) {
      force->pair->compute(eflag,vflag);
    }
    if (n_pre_reverse) {
      modify->pre_reverse(eflag,vflag);
    }
    if (force->newton) {
      comm->reverse_comm();
    }
    if (n_post_force_any) modify->post_force(vflag);
    modify->final_integrate();
    if (n_end_of_step) modify->end_of_step();
  }
}
void Verlet::cleanup()
{
  modify->post_run();
  update->update_time();
}
void Verlet::force_clear()
{
  size_t nbytes;
  if (external_force_clear) return;
  int nlocal = atom->nlocal;
  if (neighbor->includegroup == 0) {
    nbytes = sizeof(double) * nlocal;
    if (force->newton) nbytes += sizeof(double) * atom->nghost;
    if (nbytes) {
      memset(&atom->f[0][0],0,3*nbytes);
      if (torqueflag) memset(&atom->torque[0][0],0,3*nbytes);
      if (extraflag) atom->avec->force_clear(0,nbytes);
    }
  }
}
