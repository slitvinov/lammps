#include "fix_read_restart.h"
#include "atom.h"
#include "memory.h"
using namespace LAMMPS_NS;
using namespace FixConst;
FixReadRestart::FixReadRestart(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), count(nullptr), extra(nullptr)
{
  nextra = utils::inumeric(FLERR, arg[3], false, lmp);
  int nfix = utils::inumeric(FLERR, arg[4], false, lmp);
  FixReadRestart::grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);
  double **atom_extra = atom->extra;
  int nlocal = atom->nlocal;
  int i, j, m;
  for (i = 0; i < nlocal; i++) {
    m = 0;
    for (j = 0; j < nfix; j++) m += static_cast<int>(atom_extra[i][m]);
    count[i] = m;
    for (j = 0; j < m; j++) extra[i][j] = atom_extra[i][j];
  }
}
FixReadRestart::~FixReadRestart()
{
  atom->delete_callback(id, Atom::GROW);
  memory->destroy(count);
  memory->destroy(extra);
}
int FixReadRestart::setmask()
{
  int mask = 0;
  return mask;
}
double FixReadRestart::memory_usage()
{
  double bytes = (double) atom->nmax * nextra * sizeof(double);
  bytes += (double) atom->nmax * sizeof(int);
  return bytes;
}
void FixReadRestart::grow_arrays(int nmax)
{
  memory->grow(count, nmax, "read_restart:count");
  memory->grow(extra, nmax, nextra, "read_restart:extra");
}
void FixReadRestart::copy_arrays(int i, int j, int )
{
  count[j] = count[i];
  for (int m = 0; m < count[i]; m++) extra[j][m] = extra[i][m];
}
int FixReadRestart::pack_exchange(int i, double *buf)
{
  buf[0] = count[i];
  for (int m = 0; m < count[i]; m++) buf[m + 1] = extra[i][m];
  return count[i] + 1;
}
int FixReadRestart::unpack_exchange(int nlocal, double *buf)
{
  count[nlocal] = static_cast<int>(buf[0]);
  for (int m = 0; m < count[nlocal]; m++) extra[nlocal][m] = buf[m + 1];
  return count[nlocal] + 1;
}
