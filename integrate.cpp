#include "integrate.h"
#include "compute.h"
#include "force.h"
#include "modify.h"
#include "pair.h"
#include "update.h"
using namespace LAMMPS_NS;
Integrate::Integrate(LAMMPS *lmp, int , char ** ) : Pointers(lmp)
{
  elist_global = elist_atom = nullptr;
  vlist_global = vlist_atom = cvlist_atom = nullptr;
  external_force_clear = 0;
}
Integrate::~Integrate()
{
  delete [] elist_global;
  delete [] elist_atom;
  delete [] vlist_global;
  delete [] vlist_atom;
  delete [] cvlist_atom;
}
void Integrate::init()
{
  update->atimestep = update->ntimestep;
  if (force->pair && force->pair->compute_flag) pair_compute_flag = 1;
  else pair_compute_flag = 0;
}
void Integrate::ev_setup()
{
  delete [] elist_global;
  delete [] elist_atom;
  delete [] vlist_global;
  delete [] vlist_atom;
  delete [] cvlist_atom;
  elist_global = elist_atom = nullptr;
  vlist_global = vlist_atom = cvlist_atom = nullptr;
  nelist_global = nelist_atom = 0;
  nvlist_global = nvlist_atom = ncvlist_atom = 0;
  if (nelist_global) elist_global = new Compute*[nelist_global];
  if (nelist_atom) elist_atom = new Compute*[nelist_atom];
  if (nvlist_global) vlist_global = new Compute*[nvlist_global];
  if (nvlist_atom) vlist_atom = new Compute*[nvlist_atom];
  if (ncvlist_atom) cvlist_atom = new Compute*[ncvlist_atom];
  nelist_global = nelist_atom = 0;
  nvlist_global = nvlist_atom = ncvlist_atom = 0;
}
void Integrate::ev_set(bigint ntimestep)
{
  int i,flag;
  int tdflag = 0;
  flag = 0;
  int eflag_global = 0;
  for (i = 0; i < nelist_global; i++)
    if (elist_global[i]->matchstep(ntimestep)) flag = 1;
  if (flag) eflag_global = ENERGY_GLOBAL;
  flag = 0;
  int eflag_atom = 0;
  for (i = 0; i < nelist_atom; i++)
    if (elist_atom[i]->matchstep(ntimestep)) flag = 1;
  if (flag || (tdflag && nelist_atom)) eflag_atom = ENERGY_ATOM;
  if (eflag_global) update->eflag_global = ntimestep;
  if (eflag_atom) update->eflag_atom = ntimestep;
  eflag = eflag_global + eflag_atom;
  flag = 0;
  int vflag_global = 0;
  for (i = 0; i < nvlist_global; i++)
    if (vlist_global[i]->matchstep(ntimestep)) flag = 1;
  if (flag) vflag_global = virial_style;
  flag = 0;
  int vflag_atom = 0;
  for (i = 0; i < nvlist_atom; i++)
    if (vlist_atom[i]->matchstep(ntimestep)) flag = 1;
  if (flag || (tdflag && nvlist_atom)) vflag_atom = VIRIAL_ATOM;
  flag = 0;
  int cvflag_atom = 0;
  for (i = 0; i < ncvlist_atom; i++)
    if (cvlist_atom[i]->matchstep(ntimestep)) flag = 1;
  if (flag || (tdflag && ncvlist_atom)) cvflag_atom = VIRIAL_CENTROID;
  if (vflag_global) update->vflag_global = ntimestep;
  if (vflag_atom || cvflag_atom) update->vflag_atom = ntimestep;
  vflag = vflag_global + vflag_atom + cvflag_atom;
}
