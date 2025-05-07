#include <set>
#include <map>
#include <vector>
#include <cstdio>
#include <mpi.h>
#include "lammps.h"
#include "pointers.h"
#include "fix.h"
#include "fix_nve.h"
#include "atom.h"
#include "force.h"
#include "update.h"
using namespace LAMMPS_NS;
using namespace FixConst;
FixNVE::FixNVE(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg) {
  dynamic_group_allow = 1;
  time_integrate = 1;
}
int FixNVE::setmask() {
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}
void FixNVE::init() {
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
void FixNVE::initial_integrate(int) {
  double dtfm;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      dtfm = dtf / mass[type[i]];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];
    }
}
void FixNVE::final_integrate() {
  double dtfm;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      dtfm = dtf / mass[type[i]];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
    }
}
