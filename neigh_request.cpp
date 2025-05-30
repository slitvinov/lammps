#include <map>
#include <set>
#include <vector>
#include <cstdio>
#include <string>
#include <mpi.h>
#include "lammps.h"
#include "pointers.h"
#include "neigh_request.h"
#include "lmptype.h"
#include "atom.h"
#include "memory.h"
#include "neighbor.h"
using namespace LAMMPS_NS;
using namespace NeighConst;
NeighRequest::NeighRequest(LAMMPS *_lmp) : Pointers(_lmp) {
  id = 0;
  pair = 1;
  fix = compute = command = neigh = 0;
  half = 1;
  full = 0;
  occasional = 0;
  newton = 0;
  ghost = 0;
  size = 0;
  history = 0;
  granonesided = 0;
  respainner = respamiddle = respaouter = 0;
  bond = 0;
  omp = 0;
  ssa = 0;
  cut = 0;
  cutoff = 0.0;
  skip = 0;
  iskip = nullptr;
  ijskip = nullptr;
  command_style = nullptr;
  skiplist = -1;
  off2on = 0;
  copy = 0;
  trim = 0;
  copylist = -1;
  halffull = 0;
  halffulllist = -1;
  unique = 0;
  index_bin = index_pair = -1;
}
NeighRequest::NeighRequest(LAMMPS *_lmp, void *ptr, int num)
    : NeighRequest(_lmp) {
  requestor = ptr;
  requestor_instance = num;
}
NeighRequest::NeighRequest(NeighRequest *old) : NeighRequest(old->lmp) {
  copy_request(old, 1);
}
NeighRequest::~NeighRequest() {
  delete[] iskip;
  memory->destroy(ijskip);
}
void NeighRequest::copy_request(NeighRequest *other, int skipflag) {
  requestor = other->requestor;
  requestor_instance = other->requestor_instance;
  id = other->id;
  pair = other->pair;
  fix = other->fix;
  compute = other->compute;
  command = other->command;
  half = other->half;
  full = other->full;
  occasional = other->occasional;
  newton = other->newton;
  ghost = other->ghost;
  size = other->size;
  history = other->history;
  granonesided = other->granonesided;
  respainner = other->respainner;
  respamiddle = other->respamiddle;
  respaouter = other->respaouter;
  bond = other->bond;
  omp = other->omp;
  ssa = other->ssa;
  cut = other->cut;
  cutoff = other->cutoff;
  iskip = nullptr;
  ijskip = nullptr;
  if (!skipflag)
    return;
  int i, j;
  int ntp1 = atom->ntypes + 1;
  skip = other->skip;
  if (other->iskip) {
    iskip = new int[ntp1];
    for (i = 1; i < ntp1; i++)
      iskip[i] = other->iskip[i];
  }
  if (other->ijskip) {
    memory->create(ijskip, ntp1, ntp1, "neigh_request:ijskip");
    for (i = 1; i < ntp1; i++)
      for (j = 1; j < ntp1; j++)
        ijskip[i][j] = other->ijskip[i][j];
  }
}
