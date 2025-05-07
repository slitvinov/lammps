#include <map>
#include <set>
#include <unordered_set>
#include <vector>
#include <cmath>
#include <cstring>
#include <vector>
#include <string>
#include <mpi.h>
#include "lmptype.h"
#include "utils.h"
#include "lammps.h"
#include "pointers.h"
#include "neighbor.h"
#include "atom.h"
#include "atom_vec.h"
#include "pointers.h"
#include "comm.h"
#include "domain.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "nbin.h"
#include "nbin_standard.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "npair.h"
#include "pair.h"
#include "npair_half_bin_atomonly_newton.h"
#include "update.h"
using namespace LAMMPS_NS;
using namespace NeighConst;
#define RQDELTA 1
#define EXDELTA 1
#define DELTA_PERATOM 64
#define BIG 1.0e20
enum { NONE, ALL, PARTIAL, TEMPLATE };
template <typename S, typename T> static S *style_creator(LAMMPS *lmp) {
  return new T(lmp);
}
Neighbor::Neighbor(LAMMPS *lmp)
    : Pointers(lmp), pairclass(nullptr), pairnames(nullptr),
      pairmasks(nullptr) {
  MPI_Comm_rank(world, &me);
  MPI_Comm_size(world, &nprocs);
  firsttime = 1;
  style = Neighbor::BIN;
  every = 1;
  delay = 0;
  dist_check = 1;
  pgsize = 100000;
  oneatom = 2000;
  binsizeflag = 0;
  build_once = 0;
  cluster_check = 0;
  ago = -1;
  cutneighmax = 0.0;
  cutneighsq = nullptr;
  cutneighghostsq = nullptr;
  cuttype = nullptr;
  cuttypesq = nullptr;
  fixchecklist = nullptr;
  nlist = 0;
  lists = nullptr;
  nbin = 0;
  neigh_bin = nullptr;
  neigh_pair = nullptr;
  npair_perpetual = 0;
  plist = nullptr;
  nrequest = maxrequest = 0;
  requests = nullptr;
  j_sorted = nullptr;
  old_nrequest = 0;
  old_requests = nullptr;
  old_style = style;
  old_triclinic = 0;
  old_pgsize = pgsize;
  old_oneatom = oneatom;
  binclass = nullptr;
  binnames = nullptr;
  binmasks = nullptr;
  maxhold = 0;
  xhold = nullptr;
  lastcall = -1;
  last_setup_bins = -1;
  includegroup = 0;
  nex_type = maxex_type = 0;
  ex1_type = ex2_type = nullptr;
  ex_type = nullptr;
  nex_group = maxex_group = 0;
  ex1_group = ex2_group = ex1_bit = ex2_bit = nullptr;
  nex_mol = maxex_mol = 0;
  ex_mol_group = ex_mol_bit = ex_mol_intra = nullptr;
  type2collection = nullptr;
  collection2cut = nullptr;
  collection = nullptr;
  cutcollectionsq = nullptr;
  custom_collection_flag = 0;
  interval_collection_flag = 0;
  nmax_collection = 0;
  overlap_topo = 0;
}
void Neighbor::init() {
  int i, j, n;
  overlap_topo = 0;
  ncalls = ndanger = 0;
  dimension = domain->dimension;
  triclinic = domain->triclinic;
  newton_pair = force->newton_pair;
  bboxlo = domain->boxlo;
  bboxhi = domain->boxhi;
  triggersq = 0.25 * skin * skin;
  boxcheck = 0;
  n = atom->ntypes;
  if (cutneighsq == nullptr) {
    memory->create(cutneighsq, n + 1, n + 1, "neigh:cutneighsq");
    memory->create(cutneighghostsq, n + 1, n + 1, "neigh:cutneighghostsq");
    cuttype = new double[n + 1];
    cuttypesq = new double[n + 1];
  }
  double cutoff, delta, cut;
  cutneighmin = BIG;
  cutneighmax = 0.0;
  for (i = 1; i <= n; i++) {
    cuttype[i] = cuttypesq[i] = 0.0;
    for (j = 1; j <= n; j++) {
      cutoff = sqrt(force->pair->cutsq[i][j]);
      if (cutoff > 0.0)
        delta = skin;
      cut = cutoff + delta;
      cutneighsq[i][j] = cut * cut;
      cuttype[i] = MAX(cuttype[i], cut);
      cuttypesq[i] = MAX(cuttypesq[i], cut * cut);
      cutneighmin = MIN(cutneighmin, cut);
      cutneighmax = MAX(cutneighmax, cut);
      cutneighghostsq[i][j] = cut * cut;
    }
  }
  cutneighmaxsq = cutneighmax * cutneighmax;
  int respa = 0;
  delete[] fixchecklist;
  fixchecklist = nullptr;
  fixchecklist = new int[modify->nfix];
  fix_check = 0;
  memory->destroy(xhold);
  maxhold = 0;
  xhold = nullptr;
  n = atom->ntypes;
  exclude = 0;
  if (firsttime)
    init_styles();
  firsttime = 0;
  int same = init_pair();
  for (i = 0; i < nbin; i++)
    neigh_bin[i]->copy_neighbor_info();
  for (i = 0; i < nlist; i++)
    if (neigh_pair[i])
      neigh_pair[i]->copy_neighbor_info();
  for (i = 0; i < nrequest; i++) {
    delete requests[i];
    requests[i] = nullptr;
  }
  nrequest = 0;
}
void Neighbor::init_styles() {
  nbclass = 0;
#define NBIN_CLASS
#define NBinStyle(key, Class, bitmasks) nbclass++;
#include "nbin_standard.h"
#undef NBinStyle
#undef NBIN_CLASS
  binclass = new BinCreator[nbclass];
  binnames = new char *[nbclass];
  binmasks = new int[nbclass];
  nbclass = 0;
#define NBIN_CLASS
#define NBinStyle(key, Class, bitmasks)                                        \
  binnames[nbclass] = (char *)#key;                                            \
  binclass[nbclass] = &style_creator<NBin, Class>;                             \
  binmasks[nbclass++] = bitmasks;
#include "nbin_standard.h"
#undef NBinStyle
#undef NBIN_CLASS
  npclass = 0;
#define NPAIR_CLASS
#define NPairStyle(key, Class, bitmasks) npclass++;
#include "npair_half_bin_atomonly_newton.h"
#undef NPairStyle
#undef NPAIR_CLASS
  pairclass = new PairCreator[npclass];
  pairnames = new char *[npclass];
  pairmasks = new int[npclass];
  npclass = 0;
#define NPAIR_CLASS
#define NPairStyle(key, Class, bitmasks)                                       \
  pairnames[npclass] = (char *)#key;                                           \
  pairclass[npclass] = &style_creator<NPair, Class>;                           \
  pairmasks[npclass++] = bitmasks;
#include "npair_half_bin_atomonly_newton.h"
#undef NPairStyle
#undef NPAIR_CLASS
}
int Neighbor::init_pair() {
  int i, j, k, m;
  int same = 1;
  if (nrequest != old_nrequest)
    same = 0;
  if (same)
    return same;
  requests_new2old();
  delete[] lists;
  delete[] neigh_bin;
  delete[] neigh_pair;
  int nrequest_original = nrequest;
  morph_unique();
  morph_skip();
  morph_granular();
  sort_requests();
  morph_halffull();
  morph_copy_trim();
  nlist = nrequest;
  lists = new NeighList *[nrequest];
  neigh_bin = new NBin *[nrequest];
  neigh_pair = new NPair *[nrequest];
  for (i = 0; i < nrequest; i++) {
    lists[i] = new NeighList(lmp);
    lists[i]->index = i;
    lists[i]->requestor = requests[i]->requestor;
    if (requests[i]->pair) {
      lists[i]->requestor_type = NeighList::PAIR;
    } else if (requests[i]->fix) {
      lists[i]->requestor_type = NeighList::FIX;
    } else if (requests[i]->compute) {
      lists[i]->requestor_type = NeighList::COMPUTE;
    }
    if (requests[i]->pair && i < nrequest_original) {
      auto pair = (Pair *)requests[i]->requestor;
      pair->init_list(requests[i]->id, lists[i]);
    }
  }
  for (i = 0; i < nrequest; i++)
    lists[i]->post_constructor(requests[i]);
  int flag;
  for (i = 0; i < nrequest; i++) {
    flag = choose_bin(requests[i]);
    lists[i]->bin_method = flag;
    flag = choose_pair(requests[i]);
    lists[i]->pair_method = flag;
  }
  nbin = 0;
  for (i = 0; i < nrequest; i++) {
    requests[i]->index_bin = -1;
    flag = lists[i]->bin_method;
    if (flag == 0)
      continue;
    if (!requests[i]->unique) {
      for (j = 0; j < nbin; j++)
        if (neigh_bin[j]->istyle == flag && neigh_bin[j]->cutoff_custom == 0.0)
          break;
      if (j < nbin) {
        requests[i]->index_bin = j;
        continue;
      }
    }
    BinCreator &bin_creator = binclass[flag - 1];
    neigh_bin[nbin] = bin_creator(lmp);
    neigh_bin[nbin]->post_constructor(requests[i]);
    neigh_bin[nbin]->istyle = flag;
    requests[i]->index_bin = nbin;
    nbin++;
  }
  for (i = 0; i < nrequest; i++) {
    requests[i]->index_pair = -1;
    flag = lists[i]->pair_method;
    if (flag == 0) {
      neigh_pair[i] = nullptr;
      continue;
    }
    PairCreator &pair_creator = pairclass[flag - 1];
    lists[i]->np = neigh_pair[i] = pair_creator(lmp);
    neigh_pair[i]->post_constructor(requests[i]);
    neigh_pair[i]->istyle = flag;
    if (lists[i]->bin_method > 0) {
      neigh_pair[i]->nb = neigh_bin[requests[i]->index_bin];
    }
    requests[i]->index_pair = i;
  }
  for (i = 0; i < nlist; i++) {
    lists[i]->setup_pages(pgsize, oneatom);
  }
  int maxatom = atom->nmax;
  for (i = 0; i < nlist; i++) {
    if (neigh_pair[i] && (!lists[i]->copy || lists[i]->trim))
      lists[i]->grow(maxatom, maxatom);
  }
  delete[] plist;
  plist = new int[nlist];
  for (i = 0; i < nlist; i++) {
    if (lists[i]->occasional == 0 && lists[i]->pair_method)
      plist[npair_perpetual++] = i;
  }
  NeighList *ptr;
  int done = 0;
  while (!done) {
    done = 1;
    for (i = 0; i < npair_perpetual; i++) {
      for (k = 0; k < 3; k++) {
        ptr = nullptr;
        if (k == 0)
          ptr = lists[plist[i]]->listcopy;
        if (k == 1)
          ptr = lists[plist[i]]->listskip;
        if (k == 2)
          ptr = lists[plist[i]]->listfull;
        if (ptr == nullptr)
          continue;
        for (m = 0; m < nrequest; m++)
          if (ptr == lists[m])
            break;
        for (j = 0; j < npair_perpetual; j++)
          if (m == plist[j])
            break;
        if (j < i)
          continue;
        int tmp = plist[i];
        plist[i] = plist[j];
        plist[j] = tmp;
        done = 0;
        break;
      }
      if (!done)
        break;
    }
  }
  return same;
}
void Neighbor::sort_requests() {
  NeighRequest *jrq;
  int i, j, jmin;
  double jcut;
  delete[] j_sorted;
  j_sorted = new int[nrequest];
  for (i = 0; i < nrequest; i++)
    j_sorted[i] = i;
  for (i = 0; i < nrequest; i++) {
    double cutoff_min = cutneighmax;
    jmin = i;
    for (j = i; j < nrequest - 1; j++) {
      jrq = requests[j_sorted[j]];
      if (jrq->cut)
        jcut = jrq->cutoff;
      else
        jcut = cutneighmax;
      if (jcut <= cutoff_min) {
        cutoff_min = jcut;
        jmin = j;
      }
    }
    int tmp = j_sorted[i];
    j_sorted[i] = j_sorted[jmin];
    j_sorted[jmin] = tmp;
  }
}
void Neighbor::morph_unique() {
  NeighRequest *irq;
  for (int i = 0; i < nrequest; i++) {
    irq = requests[i];
    if (irq->cut) {
      if (!irq->occasional)
        irq->cutoff += skin;
      if (irq->cutoff != cutneighmax) {
        irq->unique = 1;
      } else {
        irq->cut = 0;
        irq->cutoff = 0.0;
      }
    }
  }
}
void Neighbor::morph_skip() {
  int i, j, inewton, jnewton;
  NeighRequest *irq, *jrq, *nrq;
  for (i = 0; i < nrequest; i++) {
    irq = requests[i];
    if (!irq->skip)
      continue;
  }
}
void Neighbor::morph_granular() {
  int i, j;
  NeighRequest *irq, *jrq;
  for (i = 0; i < nrequest; i++) {
    irq = requests[i];
    if (!irq->neigh)
      continue;
  }
}
void Neighbor::morph_halffull() {
  int i, j, jj;
  NeighRequest *irq, *jrq;
  double icut, jcut;
  for (i = 0; i < nrequest; i++) {
    irq = requests[i];
    int trim_flag = irq->trim;
    if (!irq->half)
      continue;
    if (irq->copy)
      continue;
    for (jj = 0; jj < nrequest; jj++) {
      if (irq->cut)
        j = j_sorted[jj];
      else
        j = jj;
      jrq = requests[j];
      if (jrq->occasional)
        continue;
      if (!jrq->full)
        continue;
      if (irq->cut)
        icut = irq->cutoff;
      else
        icut = cutneighmax;
      if (jrq->cut)
        jcut = jrq->cutoff;
      else
        jcut = cutneighmax;
      if (icut > jcut)
        continue;
      else if (icut != jcut)
        trim_flag = 1;
      if (irq->ghost != jrq->ghost)
        continue;
      if (irq->size != jrq->size)
        continue;
      if (irq->history != jrq->history)
        continue;
      if (irq->bond != jrq->bond)
        continue;
      if (irq->omp != jrq->omp)
        continue;
      if (irq->ssa != jrq->ssa)
        continue;
      if (irq->skip != jrq->skip)
        continue;
      if (irq->skip && irq->same_skip(jrq) == 0)
        continue;
      break;
    }
  }
}
void Neighbor::morph_copy_trim() {
  int i, j, jj, inewton, jnewton;
  NeighRequest *irq, *jrq;
  double icut, jcut;
  for (i = 0; i < nrequest; i++) {
    irq = requests[i];
    int trim_flag = irq->trim;
    if (irq->copy)
      continue;
    for (jj = 0; jj < nrequest; jj++) {
      if (irq->cut)
        j = j_sorted[jj];
      else
        j = jj;
      if (i == j)
        continue;
    }
  }
}
void Neighbor::requests_new2old() {
  for (int i = 0; i < old_nrequest; i++)
    delete old_requests[i];
  memory->sfree(old_requests);
  old_nrequest = nrequest;
  old_requests = (NeighRequest **)memory->smalloc(
      old_nrequest * sizeof(NeighRequest *), "neighbor:old_requests");
  for (int i = 0; i < old_nrequest; i++)
    old_requests[i] = new NeighRequest(requests[i]);
  old_style = style;
  old_triclinic = triclinic;
  old_pgsize = pgsize;
  old_oneatom = oneatom;
}
int Neighbor::choose_bin(NeighRequest *rq) {
  if (rq->skip || rq->copy || rq->halffull)
    return 0;
  int mask;
  for (int i = 0; i < nbclass; i++) {
    mask = binmasks[i];
    if (!rq->ssa != !(mask & NB_SSA))
      continue;
    if (style == Neighbor::MULTI) {
      if (!(mask & NB_MULTI))
        continue;
    } else {
      if (!(mask & NB_STANDARD))
        continue;
    }
    return i + 1;
  }
  return -1;
}
int Neighbor::choose_pair(NeighRequest *rq) {
  bool newtflag;
  newtflag = true;
  int mask;
  for (int i = 0; i < npclass; i++) {
    mask = pairmasks[i];
    return i + 1;
  }
  return -1;
}
int Neighbor::request(void *requestor, int instance) {
  if (nrequest == maxrequest) {
    maxrequest += RQDELTA;
    requests = (NeighRequest **)memory->srealloc(
        requests, maxrequest * sizeof(NeighRequest *), "neighbor:requests");
  }
  requests[nrequest] = new NeighRequest(lmp, requestor, instance);
  nrequest++;
  return nrequest - 1;
}
NeighRequest *Neighbor::add_request(Pair *requestor, int flags) {
  int irequest = request(requestor, requestor->instance_me);
  auto req = requests[irequest];
  req->apply_flags(flags);
  return req;
}
void Neighbor::setup_bins() {
  for (int i = 0; i < nbin; i++)
    neigh_bin[i]->setup_bins(style);
  last_setup_bins = update->ntimestep;
}
int Neighbor::decide() {
  ago++;
  return 1;
}
void Neighbor::build(int topoflag) {
  int i, m;
  ago = 0;
  ncalls++;
  lastcall = update->ntimestep;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  if (last_setup_bins < 0)
    setup_bins();
  for (i = 0; i < nbin; i++) {
    neigh_bin[i]->bin_atoms_setup(nall);
    neigh_bin[i]->bin_atoms();
  }
  for (i = 0; i < npair_perpetual; i++) {
    m = plist[i];
    if (!lists[m]->copy || lists[m]->trim)
      lists[m]->grow(nlocal, nall);
    neigh_pair[m]->build_setup();
    neigh_pair[m]->build(lists[m]);
  }
}
void Neighbor::set(int narg, char **arg) {
  skin = utils::numeric(FLERR, arg[0], false, lmp);
  style = Neighbor::BIN;
}
void Neighbor::modify_params(int narg, char **arg) {
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "every") == 0) {
      every = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "delay") == 0) {
      delay = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "check") == 0) {
      dist_check = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "binsize") == 0) {
      binsize_user = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      binsizeflag = 1;
      iarg += 2;
    }
  }
}
