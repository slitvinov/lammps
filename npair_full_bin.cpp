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

#include "npair_full_bin.h"

#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "error.h"
#include "my_page.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairFullBin::NPairFullBin(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   binned neighbor list construction for all neighbors
   every neighbor pair appears in list of both atoms i and j
------------------------------------------------------------------------- */

void NPairFullBin::build(NeighList *list)
{
  int i, j, k, n, itype, jtype, ibin, which, iatom;
  tagint tagprev;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq;
  int *neighptr;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;
  if (includegroup) nlocal = atom->nfirst;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;

  int inum = 0;
  ipage->reset();

  for (i = 0; i < nlocal; i++) {
    n = 0;
    neighptr = ipage->vget();

    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // loop over all atoms in surrounding bins in stencil including self
    // skip i = j

    ibin = atom2bin[i];

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status()) error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");
  }

  list->inum = inum;
  list->gnum = 0;
}
