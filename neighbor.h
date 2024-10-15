#ifndef LMP_NEIGHBOR_H
#define LMP_NEIGHBOR_H
#include "pointers.h"
namespace LAMMPS_NS {
class NeighRequest;
class NeighList;
class Neighbor : protected Pointers {
public:
  enum { NSQ, BIN, MULTI_OLD, MULTI };
  int style;
  int every;
  int delay;
  int dist_check;
  int ago;
  int pgsize;
  int oneatom;
  int includegroup;
  int build_once;
  double skin;
  double cutneighmin;
  double cutneighmax;
  double cutneighmaxsq;
  double **cutneighsq;
  double **cutneighghostsq;
  double *cuttype;
  double *cuttypesq;
  double cut_inner_sq;
  double cut_middle_sq;
  double cut_middle_inside_sq;
  int binsizeflag;
  double binsize_user;
  bigint ncalls;
  bigint ndanger;
  bigint lastcall;
  double *bboxlo, *bboxhi;
  int exclude;
  int nex_type;
  int *ex1_type, *ex2_type;
  int **ex_type;
  int nex_group;
  int *ex1_group, *ex2_group;
  int *ex1_bit, *ex2_bit;
  int nex_mol;
  int *ex_mol_group;
  int *ex_mol_bit;
  int *ex_mol_intra;
  int special_flag[4];
  int cluster_check;
  int nlist;
  int nrequest;
  int old_nrequest;
  NeighList **lists;
  NeighRequest **requests;
  NeighRequest **old_requests;
  int *j_sorted;
  int custom_collection_flag;
  int interval_collection_flag;
  int finite_cut_flag;
  int ncollections;
  int nmax_collection;
  int *type2collection;
  double *collection2cut;
  double **cutcollectionsq;
  int *collection;
  Neighbor(class LAMMPS *);
  ~Neighbor() override;
  virtual void init();
  int request(void *, int instance = 0);
  NeighRequest *add_request(class Pair *, int flags = 0);
  NeighRequest *add_request(class Fix *, int flags = 0);
  NeighRequest *add_request(class Compute *, int flags = 0);
  NeighRequest *add_request(class Command *, const char *, int flags = 0);
  int decide();
  virtual int check_distance();
  void setup_bins();
  virtual void build(int);
  void build_one(class NeighList *list, int preflag = 0);
  void set(int, char **);
  void modify_params(int, char **);
  int overlap_topo;
  NeighList *find_list(void *, const int id = 0) const;
  NeighRequest *find_request(void *, const int id = 0) const;
  const std::vector<NeighRequest *> get_pair_requests() const;
  bigint last_setup_bins;

protected:
  int me, nprocs;
  int firsttime;
  int dimension;
  int triclinic;
  int newton_pair;
  int fix_check;
  int *fixchecklist;
  double triggersq;
  double **xhold;
  int maxhold;
  int boxcheck;
  double boxlo_hold[3], boxhi_hold[3];
  double corners_hold[8][3];
  double (*corners)[3];
  double inner[2], middle[2];
  int old_style, old_triclinic;
  int old_pgsize, old_oneatom;
  int npair_perpetual;
  int *plist;
  int maxex_type;
  int maxex_group;
  int maxex_mol;
  int maxatom;
  int maxrequest;
  int nbin;
  int nbclass, npclass;
  typedef class NBin *(*BinCreator)(class LAMMPS *);
  BinCreator *binclass;
  char **binnames;
  int *binmasks;
  class NBin **neigh_bin;
  typedef class NPair *(*PairCreator)(class LAMMPS *);
  PairCreator *pairclass;
  char **pairnames;
  int *pairmasks;
  class NPair **neigh_pair;
  void init_styles();
  int init_pair();
  void sort_requests();
  void morph_unique();
  void morph_skip();
  void morph_granular();
  void morph_halffull();
  void morph_copy_trim();
  void requests_new2old();
  int choose_bin(class NeighRequest *);
  int choose_pair(class NeighRequest *);
  int copymode;
};
namespace NeighConst {
enum {
  NB_INTEL = 1 << 0,
  NB_KOKKOS_DEVICE = 1 << 1,
  NB_KOKKOS_HOST = 1 << 2,
  NB_SSA = 1 << 3,
  NB_STANDARD = 1 << 4,
  NB_MULTI = 1 << 5
};
enum {
  NS_BIN = 1 << 0,
  NS_MULTI = 1 << 1,
  NS_HALF = 1 << 2,
  NS_FULL = 1 << 3,
  NS_2D = 1 << 4,
  NS_3D = 1 << 5,
  NS_ORTHO = 1 << 6,
  NS_TRI = 1 << 7,
  NS_GHOST = 1 << 8,
  NS_SSA = 1 << 9,
  NS_MULTI_OLD = 1 << 10
};
enum {
  NP_NSQ = 1 << 0,
  NP_BIN = 1 << 1,
  NP_MULTI = 1 << 2,
  NP_HALF = 1 << 3,
  NP_FULL = 1 << 4,
  NP_ORTHO = 1 << 5,
  NP_TRI = 1 << 6,
  NP_ATOMONLY = 1 << 7,
  NP_MOLONLY = 1 << 8,
  NP_NEWTON = 1 << 9,
  NP_NEWTOFF = 1 << 10,
  NP_GHOST = 1 << 11,
  NP_SIZE = 1 << 12,
  NP_ONESIDE = 1 << 13,
  NP_RESPA = 1 << 14,
  NP_BOND = 1 << 15,
  NP_OMP = 1 << 16,
  NP_INTEL = 1 << 17,
  NP_KOKKOS_DEVICE = 1 << 18,
  NP_KOKKOS_HOST = 1 << 19,
  NP_SSA = 1 << 20,
  NP_COPY = 1 << 21,
  NP_SKIP = 1 << 22,
  NP_HALF_FULL = 1 << 23,
  NP_OFF2ON = 1 << 24,
  NP_MULTI_OLD = 1 << 25,
  NP_TRIM = 1 << 26
};
enum {
  REQ_DEFAULT = 0,
  REQ_FULL = 1 << 0,
  REQ_GHOST = 1 << 1,
  REQ_SIZE = 1 << 2,
  REQ_HISTORY = 1 << 3,
  REQ_OCCASIONAL = 1 << 4,
  REQ_RESPA_INOUT = 1 << 5,
  REQ_RESPA_ALL = 1 << 6,
  REQ_NEWTON_ON = 1 << 8,
  REQ_NEWTON_OFF = 1 << 9,
  REQ_SSA = 1 << 10,
};
} // namespace NeighConst
} // namespace LAMMPS_NS
#endif
