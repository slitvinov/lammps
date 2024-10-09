#include "fix_sdf_bounceback.h"
#include "atom.h"
#include "atom_masks.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "region.h"
#include "update.h"
#include "variable.h"
#include <stdint.h>

LAMMPS_NS::FixSdfBounceback::FixSdfBounceback(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg)
{
  arg += 3;
  narg -= 3;
  if (narg == 0) error->one(FLERR, "needs sdf file argument");
  char *sdf_path = *arg++;
  narg--;

  FILE *sdf_file;
  if ((sdf_file = fopen(sdf_path, "r")) == NULL) error->one(FLERR, "fail to open '{}'", sdf_path);

  if (fscanf(sdf_file, "%*g %*g %*g") != 0) error->one(FLERR, "fail to read '{}'", sdf_path);
  if (fscanf(sdf_file, "%" SCNd64 " %" SCNd64 " %" SCNd64 "\n", &nx_, &ny_, &nz_) != 3)
    error->one(FLERR, "fail to read '{}'", sdf_path);

  const int64_t size = nx_ * ny_ * nz_;

  sdf = (float *) memory->smalloc(size * sizeof *sdf, "fix/sdf_bounceback:sdf");

  if (fread(sdf, sizeof *sdf, size, sdf_file) != size)
    error->one(FLERR, "fail to read '{}'", sdf_path);

  if (fclose(sdf_file) != 0) error->one(FLERR, "fail to close '{}'", sdf_path);
}

LAMMPS_NS::FixSdfBounceback::~FixSdfBounceback()
{
  memory->sfree(sdf);
}

int LAMMPS_NS::FixSdfBounceback::setmask()
{
  int mask = 0;
  mask |= FixConst::POST_INTEGRATE;
  mask |= FixConst::POST_INTEGRATE_RESPA;
  return mask;
}

void LAMMPS_NS::FixSdfBounceback::post_integrate()
{
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;
  const double dt = update->dt;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      const double x1 = x[i][0];
      const double y1 = x[i][1];
      const double z1 = x[i][2];

      const double s1 = getSdf(x1, y1, z1);

      // particle is inside the domain, no need to bounce back.
      if (s1 <= 0.0) continue;

      const double x0 = x1 - dt * v[i][0];
      const double y0 = y1 - dt * v[i][1];
      const double z0 = z1 - dt * v[i][2];

      const double s0 = getSdf(x0, y0, z0);

      // particle was already outside the domain, abandon it. TODO: rescue.
      if (s0 > 0.0) continue;

      // set the particle to its previous position with reverted velocity.
      // should be placed at wall location but we use a simpler method here.
      x[i][0] = x0;
      x[i][1] = y0;
      x[i][2] = z0;

      v[i][0] *= -1;
      v[i][1] *= -1;
      v[i][2] *= -1;
    }
  }
}

float LAMMPS_NS::FixSdfBounceback::getSdf(double x, double y, double z) const
{
  // trilinear interpolation
  const double hx = domain->prd[0] / nx_;
  const double hy = domain->prd[1] / ny_;
  const double hz = domain->prd[2] / nz_;

  x -= domain->boxlo[0];
  y -= domain->boxlo[1];
  z -= domain->boxlo[2];

  int64_t ix0 = x / hx;
  int64_t iy0 = y / hy;
  int64_t iz0 = z / hz;

  const double lx = x - ix0 * hx;
  const double ly = y - iy0 * hy;
  const double lz = z - iz0 * hz;

  ix0 = (ix0 + nx_) % nx_;
  iy0 = (iy0 + ny_) % ny_;
  iz0 = (iz0 + nz_) % nz_;
  const int64_t ix1 = (ix0 + 1) % nx_;
  const int64_t iy1 = (iy0 + 1) % ny_;
  const int64_t iz1 = (iz0 + 1) % nz_;

  // if (ix0 < 0 || ix0 >= nx_) error->one(FLERR, "ix0 out of bounds '{}'", ix0);
  // if (iy0 < 0 || iy0 >= ny_) error->one(FLERR, "iy0 out of bounds '{}'", iy0);
  // if (iz0 < 0 || iz0 >= nz_) error->one(FLERR, "iz0 out of bounds '{}'", iz0);

  const double s000 = sdf[ix0 + nx_ * (iy0 + ny_ * iz0)];
  const double s100 = sdf[ix1 + nx_ * (iy0 + ny_ * iz0)];
  const double s010 = sdf[ix0 + nx_ * (iy1 + ny_ * iz0)];
  const double s110 = sdf[ix1 + nx_ * (iy1 + ny_ * iz0)];

  const double s001 = sdf[ix0 + nx_ * (iy0 + ny_ * iz1)];
  const double s101 = sdf[ix1 + nx_ * (iy0 + ny_ * iz1)];
  const double s011 = sdf[ix0 + nx_ * (iy1 + ny_ * iz1)];
  const double s111 = sdf[ix1 + nx_ * (iy1 + ny_ * iz1)];

  const double s_00 = (1.0 - lx) * s000 + lx * s100;
  const double s_10 = (1.0 - lx) * s010 + lx * s110;
  const double s_01 = (1.0 - lx) * s001 + lx * s101;
  const double s_11 = (1.0 - lx) * s011 + lx * s111;

  const double s__0 = (1.0 - ly) * s_00 + ly * s_10;
  const double s__1 = (1.0 - ly) * s_01 + ly * s_11;

  const double s___ = (1.0 - lz) * s__0 + lz * s__1;

  return s___;
}
