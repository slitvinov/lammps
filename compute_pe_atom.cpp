#include "compute_pe_atom.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "pair.h"
#include "update.h"
#include <cstring>
using namespace LAMMPS_NS;
ComputePEAtom::ComputePEAtom(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), energy(nullptr)
{
  if (narg < 3) error->all(FLERR, "Illegal compute pe/atom command");
  peratom_flag = 1;
  size_peratom_cols = 0;
  peatomflag = 1;
  timeflag = 1;
  comm_reverse = 1;
  if (narg == 3) {
    pairflag = 1;
    fixflag = 1;
  } else {
    pairflag = 0;
    fixflag = 0;
    int iarg = 3;
    while (iarg < narg) {
      if (strcmp(arg[iarg], "pair") == 0)
        pairflag = 1;
      else if (strcmp(arg[iarg], "fix") == 0)
        fixflag = 1;
      else
        error->all(FLERR, "Illegal compute pe/atom command");
      iarg++;
    }
  }
  nmax = 0;
}
ComputePEAtom::~ComputePEAtom()
{
  memory->destroy(energy);
}
void ComputePEAtom::compute_peratom()
{
  int i;
  invoked_peratom = update->ntimestep;
  if (update->eflag_atom != invoked_peratom)
    error->all(FLERR, "Per-atom energy was not tallied on needed timestep");
  if (atom->nmax > nmax) {
    memory->destroy(energy);
    nmax = atom->nmax;
    memory->create(energy, nmax, "pe/atom:energy");
    vector_atom = energy;
  }
  int nlocal = atom->nlocal;
  int npair = nlocal;
  int nbond = nlocal;
  int ntotal = nlocal;
  int nkspace = nlocal;
  if (force->newton) npair += atom->nghost;
  if (force->newton) ntotal += atom->nghost;
  for (i = 0; i < ntotal; i++) energy[i] = 0.0;
  if (pairflag && force->pair && force->pair->compute_flag) {
    double *eatom = force->pair->eatom;
    for (i = 0; i < npair; i++) energy[i] += eatom[i];
  }
  if (fixflag && modify->n_energy_atom) modify->energy_atom(nlocal, energy);
  if (force->newton) comm->reverse_comm(this);
  int *mask = atom->mask;
  for (i = 0; i < nlocal; i++)
    if (!(mask[i] & groupbit)) energy[i] = 0.0;
}
int ComputePEAtom::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m, last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = energy[i];
  return m;
}
void ComputePEAtom::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    energy[j] += buf[m++];
  }
}
double ComputePEAtom::memory_usage()
{
  double bytes = (double) nmax * sizeof(double);
  return bytes;
}
