#include "pair_dpd.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "random_mars.h"
#include "update.h"
#include <cmath>
using namespace LAMMPS_NS;
#define EPSILON 1.0e-10
PairDPD::PairDPD(LAMMPS *lmp) : Pair(lmp) {
  writedata = 1;
  random = nullptr;
}
PairDPD::~PairDPD() {
  if (copymode)
    return;
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(a0);
    memory->destroy(gamma);
    memory->destroy(sigma);
  }
  if (random)
    delete random;
}
void PairDPD::compute(int eflag, int vflag) {
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double vxtmp, vytmp, vztmp, delvx, delvy, delvz;
  double rsq, r, rinv, dot, wd, randnum, factor_dpd, factor_sqrt;
  int *ilist, *jlist, *numneigh, **firstneigh;
  evdwl = 0.0;
  ev_init(eflag, vflag);
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  double dtinvsqrt = 1.0 / sqrt(update->dt);
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    vxtmp = v[i][0];
    vytmp = v[i][1];
    vztmp = v[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_dpd = 0.0;
      factor_sqrt = special_sqrt[sbmask(j)];
      j &= NEIGHMASK;
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];
      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        if (r < EPSILON)
          continue;
        rinv = 1.0 / r;
        delvx = vxtmp - v[j][0];
        delvy = vytmp - v[j][1];
        delvz = vztmp - v[j][2];
        dot = delx * delvx + dely * delvy + delz * delvz;
        wd = 1.0 - r / cut[itype][jtype];
        randnum = random->gaussian();
        fpair = a0[itype][jtype] * wd;
        fpair -= gamma[itype][jtype] * wd * wd * dot * rinv;
        fpair *= factor_dpd;
        fpair += factor_sqrt * sigma[itype][jtype] * wd * randnum * dtinvsqrt;
        fpair *= rinv;
        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
        }
        if (eflag) {
          evdwl = 0.5 * a0[itype][jtype] * cut[itype][jtype] * wd * wd;
          evdwl *= factor_dpd;
        }
        if (evflag)
          ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely,
                   delz);
      }
    }
  }
  if (vflag_fdotr)
    virial_fdotr_compute();
}
void PairDPD::allocate() {
  int i, j;
  allocated = 1;
  int n = atom->ntypes;
  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++)
      setflag[i][j] = 0;
  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(a0, n + 1, n + 1, "pair:a0");
  memory->create(gamma, n + 1, n + 1, "pair:gamma");
  memory->create(sigma, n + 1, n + 1, "pair:sigma");
  for (i = 0; i <= atom->ntypes; i++)
    for (j = 0; j <= atom->ntypes; j++)
      sigma[i][j] = gamma[i][j] = 0.0;
}
void PairDPD::settings(int narg, char **arg) {
  if (narg != 3)
    error->all(FLERR, "Illegal pair_style command");
  temperature = utils::numeric(FLERR, arg[0], false, lmp);
  cut_global = utils::numeric(FLERR, arg[1], false, lmp);
  seed = utils::inumeric(FLERR, arg[2], false, lmp);
  if (seed <= 0)
    error->all(FLERR, "Illegal pair_style command");
  delete random;
  random = new RanMars(lmp, seed + comm->me);
  if (allocated) {
    int i, j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j])
          cut[i][j] = cut_global;
  }
}
void PairDPD::coeff(int narg, char **arg) {
  if (narg < 4 || narg > 5)
    error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated)
    allocate();
  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);
  double a0_one = utils::numeric(FLERR, arg[2], false, lmp);
  double gamma_one = utils::numeric(FLERR, arg[3], false, lmp);
  double cut_one = cut_global;
  if (narg == 5)
    cut_one = utils::numeric(FLERR, arg[4], false, lmp);
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      a0[i][j] = a0_one;
      gamma[i][j] = gamma_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }
  if (count == 0)
    error->all(FLERR, "Incorrect args for pair coefficients");
}
void PairDPD::init_style() {
  if (comm->ghost_velocity == 0)
    error->all(FLERR, "Pair dpd requires ghost atoms store velocity");
  if (force->newton_pair == 0 && comm->me == 0)
    error->warning(FLERR,
                   "Pair dpd needs newton pair on for momentum conservation");
  neighbor->add_request(this);
}
double PairDPD::init_one(int i, int j) {
  if (setflag[i][j] == 0)
    error->all(FLERR, "All pair coeffs are not set");
  sigma[i][j] = sqrt(2.0 * force->boltz * temperature * gamma[i][j]);
  cut[j][i] = cut[i][j];
  a0[j][i] = a0[i][j];
  gamma[j][i] = gamma[i][j];
  sigma[j][i] = sigma[i][j];
  return cut[i][j];
}
