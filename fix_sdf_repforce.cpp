#include "fix_sdf_repforce.h"
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

LAMMPS_NS::FixSdfRepForce::FixSdfRepForce(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  arg += 3;
  narg -= 3;
  if (narg == 0) error->one(FLERR, "needs sdf file argument");
  char *sdf_path = *arg++;
  narg--;
  if (narg == 0) error->all(FLERR, "needs scale argument");
  scale = utils::numeric(FLERR, *arg++, false, lmp);
  narg--;

  FILE *sdf_file;
  if ((sdf_file = fopen(sdf_path, "r")) == NULL) error->one(FLERR, "fail to open '{}'", sdf_path);

  if (fscanf(sdf_file, "%*g %*g %*g") != 0) error->one(FLERR, "fail to read '{}'", sdf_path);
  if (fscanf(sdf_file, "%" SCNd64 " %" SCNd64 " %" SCNd64 "\n", &nx_, &ny_, &nz_) != 3)
    error->one(FLERR, "fail to read '{}'", sdf_path);

  const int64_t size = nx_ * ny_ * nz_;

  sdf = (float *) memory->smalloc(size * sizeof *sdf, "fix/sdf_repforce:sdf");

  if (fread(sdf, sizeof *sdf, size, sdf_file) != size)
    error->one(FLERR, "fail to read '{}'", sdf_path);

  if (fclose(sdf_file) != 0) error->one(FLERR, "fail to close '{}'", sdf_path);
}

LAMMPS_NS::FixSdfRepForce::~FixSdfRepForce()
{
  memory->sfree(sdf);
}

int LAMMPS_NS::FixSdfRepForce::setmask()
{
  datamask_read = datamask_modify = 0;
  return LAMMPS_NS::FixConst::POST_FORCE;
}

void LAMMPS_NS::FixSdfRepForce::post_force(int)
{
  double **r = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      const double x = r[i][0];
      const double y = r[i][1];
      const double z = r[i][2];

      const double s = getSdf(x, y, z);

      // particle is inside the domain, force is zero
      if (s <= 0.0) continue;

      constexpr double eps = 1e-3;
      const double gx = (getSdf(x + eps, y, z) - s) / eps;
      const double gy = (getSdf(x, y + eps, z) - s) / eps;
      const double gz = (getSdf(x, y, z + eps) - s) / eps;

      const double coeff = -s * scale;

      f[i][0] += coeff * gx;
      f[i][1] += coeff * gy;
      f[i][2] += coeff * gz;
    }
  }
}

float LAMMPS_NS::FixSdfRepForce::getSdf(double x, double y, double z) const
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
