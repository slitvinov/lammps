#include <map>
#include <set>
#include <vector>
#include <cstdio>
#include <mpi.h>
#include "lammps.h"
#include "pointers.h"
#include "neigh_list.h"
#include "atom.h"
#include "comm.h"
#include "memory.h"
#include "my_page.h"
#include "neigh_request.h"
#include "neighbor.h"
using namespace LAMMPS_NS;
#define PGDELTA 1
NeighList::NeighList(LAMMPS *lmp) : Pointers(lmp) {
  maxatom = 0;
  inum = gnum = 0;
  ilist = nullptr;
  numneigh = nullptr;
  firstneigh = nullptr;
  occasional = 0;
  ghost = 0;
  ssa = 0;
  history = 0;
  respaouter = 0;
  respamiddle = 0;
  respainner = 0;
  copy = 0;
  trim = 0;
  copymode = 0;
  iskip = nullptr;
  ijskip = nullptr;
  listcopy = nullptr;
  listskip = nullptr;
  listfull = nullptr;
  fix_bond = nullptr;
  ipage = nullptr;
  inum_inner = gnum_inner = 0;
  ilist_inner = nullptr;
  numneigh_inner = nullptr;
  firstneigh_inner = nullptr;
  inum_middle = gnum_middle = 0;
  ilist_middle = nullptr;
  numneigh_middle = nullptr;
  firstneigh_middle = nullptr;
  ipage_inner = nullptr;
  ipage_middle = nullptr;
  np = nullptr;
  requestor = nullptr;
  requestor_type = NeighList::NONE;
}
NeighList::~NeighList() {
  if (copymode)
    return;
  if (!copy || trim) {
    memory->destroy(ilist);
    memory->destroy(numneigh);
    memory->sfree(firstneigh);
    delete[] ipage;
  }
  if (respainner) {
    memory->destroy(ilist_inner);
    memory->destroy(numneigh_inner);
    memory->sfree(firstneigh_inner);
    delete[] ipage_inner;
  }
  if (respamiddle) {
    memory->destroy(ilist_middle);
    memory->destroy(numneigh_middle);
    memory->sfree(firstneigh_middle);
    delete[] ipage_middle;
  }
  delete[] iskip;
  memory->destroy(ijskip);
}
void NeighList::post_constructor(NeighRequest *nq) {
  occasional = nq->occasional;
  ghost = nq->ghost;
  ssa = nq->ssa;
  history = nq->history;
  respaouter = nq->respaouter;
  respamiddle = nq->respamiddle;
  respainner = nq->respainner;
  copy = nq->copy;
  trim = nq->trim;
  id = nq->id;
  if (nq->copy) {
    listcopy = neighbor->lists[nq->copylist];
  }
  if (nq->skip) {
    listskip = neighbor->lists[nq->skiplist];
    int ntypes = atom->ntypes;
    iskip = new int[ntypes + 1];
    memory->create(ijskip, ntypes + 1, ntypes + 1, "neigh_list:ijskip");
    int i, j;
    for (i = 1; i <= ntypes; i++)
      iskip[i] = nq->iskip[i];
    for (i = 1; i <= ntypes; i++)
      for (j = 1; j <= ntypes; j++)
        ijskip[i][j] = nq->ijskip[i][j];
  }
  if (nq->halffull)
    listfull = neighbor->lists[nq->halffulllist];
  if (nq->bond)
    fix_bond = (Fix *)nq->requestor;
}
void NeighList::setup_pages(int pgsize_caller, int oneatom_caller) {
  pgsize = pgsize_caller;
  oneatom = oneatom_caller;
  int nmypage = comm->nthreads;
  ipage = new MyPage<int>[nmypage];
  for (int i = 0; i < nmypage; i++)
    ipage[i].init(oneatom, pgsize, PGDELTA);
  if (respainner) {
    ipage_inner = new MyPage<int>[nmypage];
    for (int i = 0; i < nmypage; i++)
      ipage_inner[i].init(oneatom, pgsize, PGDELTA);
  }
  if (respamiddle) {
    ipage_middle = new MyPage<int>[nmypage];
    for (int i = 0; i < nmypage; i++)
      ipage_middle[i].init(oneatom, pgsize, PGDELTA);
  }
}
void NeighList::grow(int nlocal, int nall) {
  if (ssa) {
    if ((nlocal * 3) + nall <= maxatom)
      return;
  } else if (ghost) {
    if (nall <= maxatom)
      return;
  } else {
    if (nlocal <= maxatom)
      return;
  }
  if (ssa)
    maxatom = (nlocal * 3) + nall;
  else
    maxatom = atom->nmax;
  memory->destroy(ilist);
  memory->destroy(numneigh);
  memory->sfree(firstneigh);
  memory->create(ilist, maxatom, "neighlist:ilist");
  memory->create(numneigh, maxatom, "neighlist:numneigh");
  firstneigh =
      (int **)memory->smalloc(maxatom * sizeof(int *), "neighlist:firstneigh");
  if (respainner) {
    memory->destroy(ilist_inner);
    memory->destroy(numneigh_inner);
    memory->sfree(firstneigh_inner);
    memory->create(ilist_inner, maxatom, "neighlist:ilist_inner");
    memory->create(numneigh_inner, maxatom, "neighlist:numneigh_inner");
    firstneigh_inner = (int **)memory->smalloc(maxatom * sizeof(int *),
                                               "neighlist:firstneigh_inner");
  }
  if (respamiddle) {
    memory->destroy(ilist_middle);
    memory->destroy(numneigh_middle);
    memory->sfree(firstneigh_middle);
    memory->create(ilist_middle, maxatom, "neighlist:ilist_middle");
    memory->create(numneigh_middle, maxatom, "neighlist:numneigh_middle");
    firstneigh_middle = (int **)memory->smalloc(maxatom * sizeof(int *),
                                                "neighlist:firstneigh_middle");
  }
}
