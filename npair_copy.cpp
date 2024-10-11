#include "npair_copy.h"
#include "neigh_list.h"
using namespace LAMMPS_NS;
NPairCopy::NPairCopy(LAMMPS *lmp) : NPair(lmp) {}
void NPairCopy::build(NeighList *list)
{
  NeighList *listcopy = list->listcopy;
  list->inum = listcopy->inum;
  list->gnum = listcopy->gnum;
  list->ilist = listcopy->ilist;
  list->numneigh = listcopy->numneigh;
  list->firstneigh = listcopy->firstneigh;
  list->ipage = listcopy->ipage;
}
