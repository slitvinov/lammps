/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_pe_atom.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "pair.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

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

/* ---------------------------------------------------------------------- */

ComputePEAtom::~ComputePEAtom()
{
  memory->destroy(energy);
}

/* ---------------------------------------------------------------------- */

void ComputePEAtom::compute_peratom()
{
  int i;

  invoked_peratom = update->ntimestep;
  if (update->eflag_atom != invoked_peratom)
    error->all(FLERR, "Per-atom energy was not tallied on needed timestep");

  // grow local energy array if necessary
  // needs to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(energy);
    nmax = atom->nmax;
    memory->create(energy, nmax, "pe/atom:energy");
    vector_atom = energy;
  }

  // npair includes ghosts if either newton flag is set
  //   b/c some bonds/dihedrals call pair::ev_tally with pairwise info
  // nbond includes ghosts if newton_bond is set
  // ntotal includes ghosts if either newton flag is set
  // KSpace includes ghosts if tip4pflag is set

  int nlocal = atom->nlocal;
  int npair = nlocal;
  int nbond = nlocal;
  int ntotal = nlocal;
  int nkspace = nlocal;
  if (force->newton) npair += atom->nghost;
  if (force->newton) ntotal += atom->nghost;

  // clear local energy array

  for (i = 0; i < ntotal; i++) energy[i] = 0.0;

  // add in per-atom contributions from each force

  if (pairflag && force->pair && force->pair->compute_flag) {
    double *eatom = force->pair->eatom;
    for (i = 0; i < npair; i++) energy[i] += eatom[i];
  }

  // add in per-atom contributions from relevant fixes
  // always only for owned atoms, not ghost

  if (fixflag && modify->n_energy_atom) modify->energy_atom(nlocal, energy);

  // communicate ghost energy between neighbor procs

  if (force->newton) comm->reverse_comm(this);

  // zero energy of atoms not in group
  // only do this after comm since ghost contributions must be included

  int *mask = atom->mask;

  for (i = 0; i < nlocal; i++)
    if (!(mask[i] & groupbit)) energy[i] = 0.0;
}

/* ---------------------------------------------------------------------- */

int ComputePEAtom::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = energy[i];
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputePEAtom::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    energy[j] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputePEAtom::memory_usage()
{
  double bytes = (double) nmax * sizeof(double);
  return bytes;
}
