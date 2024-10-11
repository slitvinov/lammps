#include "fix_fileforce.h"
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
struct LAMMPS_NS::FixFileforce::Force {
  int64_t n[3];
  double scale;
  float *d;
};
LAMMPS_NS::FixFileforce::FixFileforce(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  int64_t size, i, j;
  char line[2048], *force_path;
  FILE *force_file;
  float force4[4];
  arg += 3;
  narg -= 3;
  if (narg == 0) error->one(FLERR, "needs force file argument");
  force_path = *arg++;
  narg--;
  force = (struct Force *) memory->smalloc(sizeof *force, "fix/fileforce:force");
  if (narg == 0) error->all(FLERR, "needs scale argument");
  force->scale = utils::numeric(FLERR, *arg++, false, lmp);
  narg--;
  if ((force_file = fopen(force_path, "r")) == NULL)
    error->one(FLERR, "fail to open '{}'", force_path);
  if (fscanf(force_file, "%*g %*g %*g") != 0) error->one(FLERR, "fail to read '{}'", force_path);
  if (fscanf(force_file, "%" SCNd64 " %" SCNd64 " %" SCNd64 "\n", &force->n[0], &force->n[1],
             &force->n[2]) != 3)
    error->one(FLERR, "fail to read '{}'", force_path);
  size = force->n[0] * force->n[1] * force->n[2];
  force->d = (float *) memory->smalloc(3 * size * sizeof *force->d, "fix/fileforce:force");
  j = 0;
  for (i = 0; i < size; i++) {
    if (fread(force4, sizeof force4, 1, force_file) != 1)
      error->one(FLERR, "fail to read '{}'", force_path);
    force->d[j++] = force4[0];
    force->d[j++] = force4[1];
    force->d[j++] = force4[2];
  }
  if (fclose(force_file) != 0) error->one(FLERR, "fail to close '{}'", force_path);
}
LAMMPS_NS::FixFileforce::~FixFileforce()
{
  memory->sfree(force->d);
  memory->sfree(force);
}
int LAMMPS_NS::FixFileforce::setmask()
{
  datamask_read = datamask_modify = 0;
  return LAMMPS_NS::FixConst::POST_FORCE;
}
void LAMMPS_NS::FixFileforce::post_force(int)
{
  double **x = atom->x;
  double **f = atom->f;
  int nlocal, *mask = atom->mask;
  int64_t ix, iy, iz, index, size;
  imageint *image = atom->image;
  nlocal = atom->nlocal;
  size = force->n[0] * force->n[1] * force->n[2];
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      ix = force->n[0] * (x[i][0] - domain->boxlo[0]) / domain->prd[0];
      iy = force->n[1] * (x[i][1] - domain->boxlo[1]) / domain->prd[1];
      iz = force->n[2] * (x[i][2] - domain->boxlo[2]) / domain->prd[2];
      index = iz * force->n[0] * force->n[1] + iy * force->n[0] + ix;
      if (index < 0) error->one(FLERR, "index < 0");
      if (index >= size) error->one(FLERR, "index >= size");
      f[i][0] += force->scale * force->d[3 * index];
      f[i][1] += force->scale * force->d[3 * index + 1];
      f[i][2] += force->scale * force->d[3 * index + 2];
    }
}
