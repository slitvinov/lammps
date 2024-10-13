#ifndef LMP_NEIGH_LIST_H
#define LMP_NEIGH_LIST_H 
#include "pointers.h"
namespace LAMMPS_NS {
class NeighList : protected Pointers {
 public:
  enum RequestorType { NONE, PAIR, FIX, COMPUTE };
  void *requestor;
  RequestorType requestor_type;
  int index;
  int bin_method;
  int pair_method;
  int occasional;
  int ghost;
  int ssa;
  int history;
  int respaouter;
  int respamiddle;
  int respainner;
  int copy;
  int trim;
  int copymode;
  int id;
  int inum;
  int gnum;
  int *ilist;
  int *numneigh;
  int **firstneigh;
  int maxatom;
  int pgsize;
  int oneatom;
  MyPage<int> *ipage;
  int inum_inner;
  int gnum_inner;
  int *ilist_inner;
  int *numneigh_inner;
  int **firstneigh_inner;
  int inum_middle;
  int gnum_middle;
  int *ilist_middle;
  int *numneigh_middle;
  int **firstneigh_middle;
  MyPage<int> *ipage_inner;
  MyPage<int> *ipage_middle;
  int *iskip;
  int **ijskip;
  NeighList *listcopy;
  NeighList *listskip;
  NeighList *listfull;
  class Fix *fix_bond;
  class NPair *np;
  NeighList(class LAMMPS *);
  ~NeighList() override;
  void post_constructor(class NeighRequest *);
  void setup_pages(int, int);
  void grow(int, int);
  void print_attributes();
  int get_maxlocal() { return maxatom; }
};
}
#endif
