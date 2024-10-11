#include "neighbor.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "nbin.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "npair.h"
#include "pair.h"
#include "style_nbin.h"
#include "style_npair.h"
#include "suffix.h"
#include "tokenizer.h"
#include "update.h"
#include <cmath>
#include <cstring>
using namespace LAMMPS_NS;
using namespace NeighConst;
#define RQDELTA 1
#define EXDELTA 1
#define DELTA_PERATOM 64
#define BIG 1.0e20
enum{NONE,ALL,PARTIAL,TEMPLATE};
static const char cite_neigh_multi_old[] =
  "neighbor multi/old command: doi:10.1016/j.cpc.2008.03.005\n\n"
  "@Article{Intveld08,\n"
  " author =  {in 't Veld, P. J. and S. J. Plimpton and G. S. Grest},\n"
  " title =   {Accurate and Efficient Methods for Modeling Colloidal\n"
  "            Mixtures in an Explicit Solvent using Molecular Dynamics},\n"
  " journal = {Comput.\\ Phys.\\ Commun.},\n"
  " year =    2008,\n"
  " volume =  179,\n"
  " number =  5,\n"
  " pages =   {320--329}\n"
  "}\n\n";
static const char cite_neigh_multi[] =
  "neighbor multi command: doi:10.1016/j.cpc.2008.03.005, doi:10.1007/s40571-020-00361-2\n\n"
  "@Article{Intveld08,\n"
  " author =  {in 't Veld, P. J. and S. J.~Plimpton and G. S. Grest},\n"
  " title =   {Accurate and Efficient Methods for Modeling Colloidal\n"
  "            Mixtures in an Explicit Solvent using Molecular Dynamics},\n"
  " journal = {Comput.\\ Phys.\\ Commut.},\n"
  " year =    2008,\n"
  " volume =  179,\n"
  " pages =   {320--329}\n"
  "}\n\n"
  "@article{Shire2020,\n"
  " author = {Shire, Tom and Hanley, Kevin J. and Stratford, Kevin},\n"
  " title = {{DEM} Simulations of Polydisperse Media: Efficient Contact\n"
  "          Detection Applied to Investigate the Quasi-Static Limit},\n"
  " journal = {Computational Particle Mechanics},\n"
  " year = {2020}\n"
  "@article{Monti2022,\n"
  " author = {Monti, Joseph M. and Clemmer, Joel T. and Srivastava, \n"
  "           Ishan and Silbert, Leonardo E. and Grest, Gary S. \n"
  "           and Lechman, Jeremy B.},\n"
  " title = {Large-scale frictionless jamming with power-law particle \n"
  "          size distributions},\n"
  " journal = {Phys. Rev. E},\n"
  " volume = {106}\n"
  " issue = {3}\n"
  " year = {2022}\n"
  "}\n\n";
template <typename S, typename T> static S *style_creator(LAMMPS *lmp)
{
  return new T(lmp);
}
Neighbor::Neighbor(LAMMPS *lmp) : Pointers(lmp),
pairclass(nullptr), pairnames(nullptr), pairmasks(nullptr)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
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
  copymode = 0;
  overlap_topo = 0;
}
Neighbor::~Neighbor()
{
  if (copymode) return;
  memory->destroy(cutneighsq);
  memory->destroy(cutneighghostsq);
  delete[] cuttype;
  delete[] cuttypesq;
  delete[] fixchecklist;
  for (int i = 0; i < nlist; i++) delete lists[i];
  for (int i = 0; i < nbin; i++) delete neigh_bin[i];
  for (int i = 0; i < nlist; i++) delete neigh_pair[i];
  delete[] lists;
  delete[] neigh_bin;
  delete[] neigh_pair;
  delete[] plist;
  for (int i = 0; i < nrequest; i++)
    if (requests[i]) delete requests[i];
  memory->sfree(requests);
  for (int i = 0; i < old_nrequest; i++)
    if (old_requests[i]) delete old_requests[i];
  memory->sfree(old_requests);
  delete[] j_sorted;
  delete[] binclass;
  delete[] binnames;
  delete[] binmasks;
  delete[] pairclass;
  delete[] pairnames;
  delete[] pairmasks;
  memory->destroy(xhold);
  memory->destroy(ex1_type);
  memory->destroy(ex2_type);
  memory->destroy(ex_type);
  memory->destroy(ex1_group);
  memory->destroy(ex2_group);
  delete[] ex1_bit;
  delete[] ex2_bit;
  memory->destroy(ex_mol_group);
  delete[] ex_mol_bit;
  memory->destroy(ex_mol_intra);
  memory->destroy(type2collection);
  memory->destroy(collection2cut);
  memory->destroy(collection);
  memory->destroy(cutcollectionsq);
}
void Neighbor::init()
{
  int i,j,n;
  overlap_topo = 0;
  ncalls = ndanger = 0;
  dimension = domain->dimension;
  triclinic = domain->triclinic;
  newton_pair = force->newton_pair;
  if (delay > 0 && (delay % every) != 0)
    error->all(FLERR,"Neighbor delay must be 0 or multiple of every setting");
  if (pgsize < 10*oneatom)
    error->all(FLERR,"Neighbor page size must be >= 10x the one atom setting");
  if (triclinic == 0) {
    bboxlo = domain->boxlo;
    bboxhi = domain->boxhi;
  } else {
    bboxlo = domain->boxlo_bound;
    bboxhi = domain->boxhi_bound;
  }
  triggersq = 0.25*skin*skin;
  boxcheck = 0;
  if (domain->box_change && (domain->xperiodic || domain->yperiodic ||
                             (dimension == 3 && domain->zperiodic)))
      boxcheck = 1;
  n = atom->ntypes;
  if (cutneighsq == nullptr) {
    memory->create(cutneighsq,n+1,n+1,"neigh:cutneighsq");
    memory->create(cutneighghostsq,n+1,n+1,"neigh:cutneighghostsq");
    cuttype = new double[n+1];
    cuttypesq = new double[n+1];
  }
  double cutoff,delta,cut;
  cutneighmin = BIG;
  cutneighmax = 0.0;
  for (i = 1; i <= n; i++) {
    cuttype[i] = cuttypesq[i] = 0.0;
    for (j = 1; j <= n; j++) {
      if (force->pair) cutoff = sqrt(force->pair->cutsq[i][j]);
      else cutoff = 0.0;
      if (cutoff > 0.0) delta = skin;
      else delta = 0.0;
      cut = cutoff + delta;
      cutneighsq[i][j] = cut*cut;
      cuttype[i] = MAX(cuttype[i],cut);
      cuttypesq[i] = MAX(cuttypesq[i],cut*cut);
      cutneighmin = MIN(cutneighmin,cut);
      cutneighmax = MAX(cutneighmax,cut);
      if (force->pair && force->pair->ghostneigh) {
        cut = force->pair->cutghost[i][j] + skin;
        cutneighghostsq[i][j] = cut*cut;
      } else cutneighghostsq[i][j] = cut*cut;
    }
  }
  cutneighmaxsq = cutneighmax * cutneighmax;
  if (style == Neighbor::MULTI) {
    int icollection, jcollection;
    if (!custom_collection_flag) {
      ncollections = n;
      interval_collection_flag = 0;
      if (!type2collection)
        memory->create(type2collection,n+1,"neigh:type2collection");
      for (i = 1; i <= n; i++)
        type2collection[i] = i-1;
    }
    memory->grow(cutcollectionsq, ncollections, ncollections, "neigh:cutcollectionsq");
    for (i = 0; i < ncollections; i++)
      for (j = 0; j < ncollections; j++)
        cutcollectionsq[i][j] = 0.0;
    if (!interval_collection_flag) {
      finite_cut_flag = 0;
      for (i = 1; i <= n; i++){
        icollection = type2collection[i];
        for (j = 1; j <= n; j++){
          jcollection = type2collection[j];
          if (cutneighsq[i][j] > cutcollectionsq[icollection][jcollection]) {
            cutcollectionsq[icollection][jcollection] = cutneighsq[i][j];
            cutcollectionsq[jcollection][icollection] = cutneighsq[i][j];
          }
        }
      }
    } else {
      if (force->pair->finitecutflag) {
        finite_cut_flag = 1;
        double ri, rj, tmp;
        for (i = 0; i < ncollections; i++){
          ri = collection2cut[i]*0.5;
          for (j = 0; j < ncollections; j++){
            rj = collection2cut[j]*0.5;
            tmp = force->pair->radii2cut(ri, rj) + skin;
            cutcollectionsq[i][j] = tmp*tmp;
          }
        }
      } else {
        finite_cut_flag = 0;
        if (!type2collection)
          memory->create(type2collection,n+1,"neigh:type2collection");
        for (i = 1; i <= n; i++)
          type2collection[i] = -1;
        double cuttmp;
        for (i = 1; i <= n; i++){
          cuttmp = sqrt(cutneighsq[i][i]) - skin;
          for (icollection = 0; icollection < ncollections; icollection ++){
            if (collection2cut[icollection] >= cuttmp) {
              type2collection[i] = icollection;
              break;
            }
          }
          if (type2collection[i] == -1)
            error->all(FLERR, "Pair cutoff exceeds interval cutoffs for multi");
        }
        for (i = 1; i <= n; i++){
          icollection = type2collection[i];
          for (j = 1; j <= n; j++){
            jcollection = type2collection[j];
            if (cutneighsq[i][j] > cutcollectionsq[icollection][jcollection]) {
              cutcollectionsq[icollection][jcollection] = cutneighsq[i][j];
              cutcollectionsq[jcollection][icollection] = cutneighsq[i][j];
            }
          }
        }
      }
    }
  }
  int respa = 0;
  delete[] fixchecklist;
  fixchecklist = nullptr;
  fixchecklist = new int[modify->nfix];
  fix_check = 0;
  for (i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->force_reneighbor)
      fixchecklist[fix_check++] = i;
  if (dist_check == 0) {
    memory->destroy(xhold);
    maxhold = 0;
    xhold = nullptr;
  }
  if (dist_check) {
    if (maxhold == 0) {
      maxhold = atom->nmax;
      memory->create(xhold,maxhold,3,"neigh:xhold");
    }
  }
  n = atom->ntypes;
  if (nex_type == 0 && nex_group == 0 && nex_mol == 0) exclude = 0;
  else exclude = 1;
  if (nex_type) {
    memory->destroy(ex_type);
    memory->create(ex_type,n+1,n+1,"neigh:ex_type");
    for (i = 1; i <= n; i++)
      for (j = 1; j <= n; j++)
        ex_type[i][j] = 0;
    for (i = 0; i < nex_type; i++) {
      if (ex1_type[i] <= 0 || ex1_type[i] > n ||
          ex2_type[i] <= 0 || ex2_type[i] > n)
        error->all(FLERR,"Invalid atom type in neighbor exclusion list");
      ex_type[ex1_type[i]][ex2_type[i]] = 1;
      ex_type[ex2_type[i]][ex1_type[i]] = 1;
    }
  }
  if (nex_group) {
    delete[] ex1_bit;
    delete[] ex2_bit;
    ex1_bit = new int[nex_group];
    ex2_bit = new int[nex_group];
    for (i = 0; i < nex_group; i++) {
      ex1_bit[i] = group->bitmask[ex1_group[i]];
      ex2_bit[i] = group->bitmask[ex2_group[i]];
    }
  }
  if (nex_mol) {
    delete[] ex_mol_bit;
    ex_mol_bit = new int[nex_mol];
    for (i = 0; i < nex_mol; i++)
      ex_mol_bit[i] = group->bitmask[ex_mol_group[i]];
  }
  if (firsttime) init_styles();
  firsttime = 0;
  int same = init_pair();
  for (i = 0; i < nbin; i++) neigh_bin[i]->copy_neighbor_info();
  for (i = 0; i < nlist; i++)
    if (neigh_pair[i]) neigh_pair[i]->copy_neighbor_info();
  if (!same && (nrequest > 0) && (comm->me == 0)) print_pairwise_info();
  for (i = 0; i < nrequest; i++) {
    delete requests[i];
    requests[i] = nullptr;
  }
  nrequest = 0;
}
void Neighbor::init_styles()
{
  nbclass = 0;
#define NBIN_CLASS 
#define NBinStyle(key,Class,bitmasks) nbclass++;
#include "style_nbin.h"
#undef NBinStyle
#undef NBIN_CLASS
  binclass = new BinCreator[nbclass];
  binnames = new char*[nbclass];
  binmasks = new int[nbclass];
  nbclass = 0;
#define NBIN_CLASS 
#define NBinStyle(key,Class,bitmasks) \
  binnames[nbclass] = (char *) #key; \
  binclass[nbclass] = &style_creator<NBin, Class>; \
  binmasks[nbclass++] = bitmasks;
#include "style_nbin.h"
#undef NBinStyle
#undef NBIN_CLASS
  npclass = 0;
#define NPAIR_CLASS 
#define NPairStyle(key,Class,bitmasks) npclass++;
#include "style_npair.h"
#undef NPairStyle
#undef NPAIR_CLASS
  pairclass = new PairCreator[npclass];
  pairnames = new char*[npclass];
  pairmasks = new int[npclass];
  npclass = 0;
#define NPAIR_CLASS 
#define NPairStyle(key,Class,bitmasks) \
  pairnames[npclass] = (char *) #key; \
  pairclass[npclass] = &style_creator<NPair, Class>; \
  pairmasks[npclass++] = bitmasks;
#include "style_npair.h"
#undef NPairStyle
#undef NPAIR_CLASS
}
int Neighbor::init_pair()
{
  int i,j,k,m;
  int same = 1;
  if (style != old_style) same = 0;
  if (triclinic != old_triclinic) same = 0;
  if (pgsize != old_pgsize) same = 0;
  if (oneatom != old_oneatom) same = 0;
  if (nrequest != old_nrequest) same = 0;
  else
    for (i = 0; i < nrequest; i++)
      if (requests[i]->identical(old_requests[i]) == 0) same = 0;
#ifdef NEIGH_LIST_DEBUG
  if (comm->me == 0) printf("SAME flag %d\n",same);
#endif
  if (same) return same;
  requests_new2old();
  for (i = 0; i < nlist; i++) delete lists[i];
  for (i = 0; i < nbin; i++) delete neigh_bin[i];
  for (i = 0; i < nlist; i++) delete neigh_pair[i];
  delete[] lists;
  delete[] neigh_bin;
  delete[] neigh_pair;
  if (style == Neighbor::BIN) {
    for (i = 0; i < nrequest; i++)
      if (requests[i]->occasional && requests[i]->ghost)
        error->all(FLERR,"Cannot request an occasional binned neighbor list "
                   "with ghost info");
  }
  int nrequest_original = nrequest;
  morph_unique();
  morph_skip();
  morph_granular();
  sort_requests();
  morph_halffull();
  morph_copy_trim();
  nlist = nrequest;
  lists = new NeighList*[nrequest];
  neigh_bin = new NBin*[nrequest];
  neigh_pair = new NPair*[nrequest];
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
      auto pair = (Pair *) requests[i]->requestor;
      pair->init_list(requests[i]->id,lists[i]);
    } else if (requests[i]->fix && i < nrequest_original) {
      Fix *fix = (Fix *) requests[i]->requestor;
      fix->init_list(requests[i]->id,lists[i]);
    } else if (requests[i]->compute && i < nrequest_original) {
      auto compute = (Compute *) requests[i]->requestor;
      compute->init_list(requests[i]->id,lists[i]);
    }
  }
  for (i = 0; i < nrequest; i++)
    lists[i]->post_constructor(requests[i]);
  int flag;
  for (i = 0; i < nrequest; i++) {
    flag = choose_bin(requests[i]);
    lists[i]->bin_method = flag;
    if (flag < 0)
      error->all(FLERR,"Requested neighbor bin option does not exist");
    flag = choose_pair(requests[i]);
    lists[i]->pair_method = flag;
    if (flag < 0)
      error->all(FLERR,"Requested neighbor pair method does not exist");
  }
  nbin = 0;
  for (i = 0; i < nrequest; i++) {
    requests[i]->index_bin = -1;
    flag = lists[i]->bin_method;
    if (flag == 0) continue;
    if (!requests[i]->unique) {
      for (j = 0; j < nbin; j++)
        if (neigh_bin[j]->istyle == flag &&
            neigh_bin[j]->cutoff_custom == 0.0) break;
      if (j < nbin) {
        requests[i]->index_bin = j;
        continue;
      }
    }
    BinCreator &bin_creator = binclass[flag-1];
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
    PairCreator &pair_creator = pairclass[flag-1];
    lists[i]->np = neigh_pair[i] = pair_creator(lmp);
    neigh_pair[i]->post_constructor(requests[i]);
    neigh_pair[i]->istyle = flag;
    if (lists[i]->bin_method > 0) {
      neigh_pair[i]->nb = neigh_bin[requests[i]->index_bin];
      if (neigh_pair[i]->nb == nullptr)
        error->all(FLERR,"Could not assign bin method to neighbor pair");
    }
    requests[i]->index_pair = i;
  }
  for (i = 0; i < nlist; i++) {
    lists[i]->setup_pages(pgsize,oneatom);
  }
  int maxatom = atom->nmax;
  for (i = 0; i < nlist; i++) {
    if (neigh_pair[i] && (!lists[i]->copy || lists[i]->trim))
      lists[i]->grow(maxatom,maxatom);
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
        if (k == 0) ptr = lists[plist[i]]->listcopy;
        if (k == 1) ptr = lists[plist[i]]->listskip;
        if (k == 2) ptr = lists[plist[i]]->listfull;
        if (ptr == nullptr) continue;
        for (m = 0; m < nrequest; m++)
          if (ptr == lists[m]) break;
        for (j = 0; j < npair_perpetual; j++)
          if (m == plist[j]) break;
        if (j < i) continue;
        int tmp = plist[i];
        plist[i] = plist[j];
        plist[j] = tmp;
        done = 0;
        break;
      }
      if (!done) break;
    }
  }
#ifdef NEIGH_LIST_DEBUG
  for (i = 0; i < nrequest; i++) lists[i]->print_attributes();
#endif
  return same;
}
void Neighbor::sort_requests()
{
  NeighRequest *jrq;
  int i,j,jmin;
  double jcut;
  delete[] j_sorted;
  j_sorted = new int[nrequest];
  for (i = 0; i < nrequest; i++)
    j_sorted[i] = i;
  for (i = 0; i < nrequest; i++) {
    double cutoff_min = cutneighmax;
    jmin = i;
    for (j = i; j < nrequest-1; j++) {
      jrq = requests[j_sorted[j]];
      if (jrq->cut) jcut = jrq->cutoff;
      else jcut = cutneighmax;
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
void Neighbor::morph_unique()
{
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
void Neighbor::morph_skip()
{
  int i,j,inewton,jnewton;
  NeighRequest *irq,*jrq,*nrq;
  for (i = 0; i < nrequest; i++) {
    irq = requests[i];
    if (!irq->skip) continue;
    if (irq->halffull) continue;
    if (irq->copy) continue;
    for (j = 0; j < nrequest; j++) {
      if (i == j) continue;
      jrq = requests[j];
      if (jrq->occasional) continue;
      if (jrq->skip) continue;
      if (irq->half != jrq->half) continue;
      if (irq->full != jrq->full) continue;
      inewton = irq->newton;
      if (inewton == 0) inewton = force->newton_pair ? 1 : 2;
      jnewton = jrq->newton;
      if (jnewton == 0) jnewton = force->newton_pair ? 1 : 2;
      if (inewton != jnewton) continue;
      if (irq->ghost != jrq->ghost) continue;
      if (irq->size != jrq->size) continue;
      if (irq->history != jrq->history) continue;
      if (irq->bond != jrq->bond) continue;
      if (irq->omp != jrq->omp) continue;
      if (irq->ssa != jrq->ssa) continue;
      if (irq->cut != jrq->cut) continue;
      if (irq->cutoff != jrq->cutoff) continue;
      break;
    }
    if (j < nrequest) irq->skiplist = j;
    else {
      int newrequest = request(this,-1);
      irq->skiplist = newrequest;
      nrq = requests[newrequest];
      nrq->copy_request(irq,0);
      nrq->pair = nrq->fix = nrq->compute = nrq->command = 0;
      nrq->neigh = 1;
      nrq->skip = 0;
      if (irq->unique) nrq->unique = 1;
    }
  }
}
void Neighbor::morph_granular()
{
  int i,j;
  NeighRequest *irq,*jrq;
  for (i = 0; i < nrequest; i++) {
    irq = requests[i];
    if (!irq->neigh) continue;
    if (!irq->size) continue;
    int onesided = -1;
    for (j = 0; j < nrequest; j++) {
      jrq = requests[j];
      if (!jrq->pair) continue;
      if (!jrq->size) continue;
      if (!jrq->skip || jrq->skiplist != i) continue;
      if (onesided < 0) onesided = jrq->granonesided;
      else if (onesided != jrq->granonesided) onesided = 2;
      if (onesided == 2) break;
    }
    if (onesided == 2) {
      irq->newton = 2;
      irq->granonesided = 0;
      for (j = 0; j < nrequest; j++) {
        jrq = requests[j];
        if (!jrq->pair) continue;
        if (!jrq->size) continue;
        if (!jrq->skip || jrq->skiplist != i) continue;
        jrq->off2on = 1;
      }
    }
  }
}
void Neighbor::morph_halffull()
{
  int i,j,jj;
  NeighRequest *irq,*jrq;
  double icut,jcut;
  for (i = 0; i < nrequest; i++) {
    irq = requests[i];
    int trim_flag = irq->trim;
    if (!irq->half) continue;
    if (irq->copy) continue;
    for (jj = 0; jj < nrequest; jj++) {
      if (irq->cut) j = j_sorted[jj];
      else j = jj;
      jrq = requests[j];
      if (jrq->occasional) continue;
      if (!jrq->full) continue;
      if (irq->cut) icut = irq->cutoff;
      else icut = cutneighmax;
      if (jrq->cut) jcut = jrq->cutoff;
      else jcut = cutneighmax;
      if (icut > jcut) continue;
      else if (icut != jcut) trim_flag = 1;
      if (irq->ghost != jrq->ghost) continue;
      if (irq->size != jrq->size) continue;
      if (irq->history != jrq->history) continue;
      if (irq->bond != jrq->bond) continue;
      if (irq->omp != jrq->omp) continue;
      if (irq->ssa != jrq->ssa) continue;
      if (irq->skip != jrq->skip) continue;
      if (irq->skip && irq->same_skip(jrq) == 0) continue;
      break;
    }
    if (jj < nrequest) {
      irq->halffull = 1;
      irq->halffulllist = j;
      irq->trim = trim_flag;
    }
  }
}
void Neighbor::morph_copy_trim()
{
  int i,j,jj,inewton,jnewton;
  NeighRequest *irq,*jrq;
  double icut,jcut;
  for (i = 0; i < nrequest; i++) {
    irq = requests[i];
    int trim_flag = irq->trim;
    if (irq->copy) continue;
    for (jj = 0; jj < nrequest; jj++) {
      if (irq->cut) j = j_sorted[jj];
      else j = jj;
      if (i == j) continue;
      jrq = requests[j];
      if (jrq->copy && jrq->copylist == i) continue;
      if (irq->cut) icut = irq->cutoff;
      else icut = cutneighmax;
      if (jrq->cut) jcut = jrq->cutoff;
      else jcut = cutneighmax;
      if (icut > jcut) continue;
      else if (icut != jcut) trim_flag = 1;
      if (jrq->occasional) continue;
      if (!irq->occasional && !irq->cut && j > i) continue;
      if (irq->half != jrq->half) continue;
      if (irq->full != jrq->full) continue;
      inewton = irq->newton;
      if (inewton == 0) inewton = force->newton_pair ? 1 : 2;
      jnewton = jrq->newton;
      if (jnewton == 0) jnewton = force->newton_pair ? 1 : 2;
      if (inewton != jnewton) continue;
      if (irq->ghost && !jrq->ghost) continue;
      if (jrq->respamiddle) continue;
      if (jrq->respainner) continue;
      if (irq->size != jrq->size) continue;
      if (irq->history != jrq->history) continue;
      if (irq->bond != jrq->bond) continue;
      if (irq->ssa != jrq->ssa) continue;
      if (irq->skip != jrq->skip) continue;
      if (irq->skip && irq->same_skip(jrq) == 0) continue;
      break;
    }
    if (jj < nrequest) {
      irq->copy = 1;
      irq->trim = trim_flag;
      if (jrq->copy && irq->cutoff == requests[jrq->copylist]->cutoff)
        irq->copylist = jrq->copylist;
      else
        irq->copylist = j;
    }
  }
}
void Neighbor::print_pairwise_info()
{
  int i;
  NeighRequest *rq;
  const double cutghost = MAX(cutneighmax,comm->cutghostuser);
  double binsize, bbox[3];
  bbox[0] = bboxhi[0]-bboxlo[0];
  bbox[1] = bboxhi[1]-bboxlo[1];
  bbox[2] = bboxhi[2]-bboxlo[2];
  if (binsizeflag) binsize = binsize_user;
  else if (style == Neighbor::BIN) binsize = 0.5*cutneighmax;
  else binsize = 0.5*cutneighmin;
  if (binsize == 0.0) binsize = bbox[0];
  int nperpetual = 0;
  int noccasional = 0;
  int nextra = 0;
  for (i = 0; i < nlist; i++) {
    if (lists[i]->pair_method == 0) nextra++;
    else if (lists[i]->occasional) noccasional++;
    else nperpetual++;
  }
  std::string out = "Neighbor list info ...\n";
  out += fmt::format("  update: every = {} steps, delay = {} steps, check = {}\n",
                     every,delay,dist_check ? "yes" : "no");
  out += fmt::format("  max neighbors/atom: {}, page size: {}\n",
                     oneatom, pgsize);
  out += fmt::format("  master list distance cutoff = {:.8g}\n",cutneighmax);
  out += fmt::format("  ghost atom cutoff = {:.8g}\n",cutghost);
  if (style != Neighbor::NSQ)
    out += fmt::format("  binsize = {:.8g}, bins = {:g} {:g} {:g}\n",binsize,
                       ceil(bbox[0]/binsize), ceil(bbox[1]/binsize),
                       ceil(bbox[2]/binsize));
  out += fmt::format("  {} neighbor lists, perpetual/occasional/extra = {} {} {}\n",
                     nlist,nperpetual,noccasional,nextra);
  for (i = 0; i < nlist; i++) {
    rq = requests[i];
    if (rq->pair) {
      char *pname = force->pair_match_ptr((Pair *) rq->requestor);
      if (pname) out += fmt::format("  ({}) pair {}",i+1,pname);
      else out += fmt::format("  ({}) pair (none)",i+1);
    } else if (rq->fix) {
      out += fmt::format("  ({}) fix {}",i+1,((Fix *) rq->requestor)->style);
    } else if (rq->compute) {
      out += fmt::format("  ({}) compute {}",i+1,((Compute *) rq->requestor)->style);
    } else if (rq->command) {
      out += fmt::format("  ({}) command {}",i+1,rq->command_style);
    } else if (rq->neigh) {
      out += fmt::format("  ({}) neighbor class addition",i+1);
    }
    if (rq->occasional) out += ", occasional";
    else out += ", perpetual";
    if (rq->copy) {
      if (rq->trim)
        out += fmt::format(", trim from ({})",rq->copylist+1);
      else
        out += fmt::format(", copy from ({})",rq->copylist+1);
    } else if (rq->halffull)
      if (rq->trim)
        out += fmt::format(", half/full trim from ({})",rq->halffulllist+1);
      else
        out += fmt::format(", half/full from ({})",rq->halffulllist+1);
    else if (rq->skip)
      out += fmt::format(", skip from ({})",rq->skiplist+1);
    out += "\n";
    out += "      attributes: ";
    if (rq->half) out += "half";
    else if (rq->full) out += "full";
    if (rq->newton == 0) {
      if (force->newton_pair) out += ", newton on";
      else out += ", newton off";
    } else if (rq->newton == 1) out += ", newton on";
    else if (rq->newton == 2) out += ", newton off";
    if (rq->ghost) out += ", ghost";
    if (rq->size) out += ", size";
    if (rq->history) out += ", history";
    if (rq->granonesided) out += ", onesided";
    if (rq->respamiddle) out += ", respa outer/middle/inner";
    else if (rq->respainner) out += ", respa outer/inner";
    if (rq->bond) out += ", bond";
    if (rq->omp) out += ", omp";
    if (rq->ssa) out += ", ssa";
    if (rq->cut) out += fmt::format(", cut {}",rq->cutoff);
    if (rq->off2on) out += ", off2on";
    out += "\n";
    out += "      ";
    if (lists[i]->pair_method == 0) out += "pair build: none\n";
    else out += fmt::format("pair build: {}\n",pairnames[lists[i]->pair_method-1]);
    out += "      ";
    if (lists[i]->bin_method == 0) out += "bin: none\n";
    else out += fmt::format("bin: {}\n",binnames[lists[i]->bin_method-1]);
  }
  utils::logmesg(lmp,out);
}
void Neighbor::requests_new2old()
{
  for (int i = 0; i < old_nrequest; i++) delete old_requests[i];
  memory->sfree(old_requests);
  old_nrequest = nrequest;
  old_requests = (NeighRequest **)
    memory->smalloc(old_nrequest*sizeof(NeighRequest *),"neighbor:old_requests");
  for (int i = 0; i < old_nrequest; i++)
    old_requests[i] = new NeighRequest(requests[i]);
  old_style = style;
  old_triclinic = triclinic;
  old_pgsize = pgsize;
  old_oneatom = oneatom;
}
NeighRequest *Neighbor::find_request(void *classptr, const int id) const
{
  if (classptr == nullptr) return nullptr;
  for (int i = 0; i < nrequest; i++)
    if ((requests[i]->requestor == classptr) && (requests[i]->id == id)) return requests[i];
  return nullptr;
}
const std::vector<NeighRequest *> Neighbor::get_pair_requests() const
{
  std::vector<NeighRequest *> matches;
  for (int i=0; i < nrequest; ++i)
    if (requests[i]->pair) matches.push_back(requests[i]);
  return matches;
}
NeighList *Neighbor::find_list(void *classptr, const int id) const
{
  if (classptr == nullptr) return nullptr;
  for (int i = 0; i < nlist; i++)
    if ((lists[i]->requestor == classptr) && (lists[i]->id == id)) return lists[i];
  return nullptr;
}
int Neighbor::choose_bin(NeighRequest *rq)
{
  if (style == Neighbor::NSQ) return 0;
  if (rq->skip || rq->copy || rq->halffull) return 0;
  int mask;
  for (int i = 0; i < nbclass; i++) {
    mask = binmasks[i];
    if (!rq->ssa != !(mask & NB_SSA)) continue;
    if (style == Neighbor::MULTI) {
      if (!(mask & NB_MULTI)) continue;
    } else {
      if (!(mask & NB_STANDARD)) continue;
    }
    return i+1;
  }
  return -1;
}
int Neighbor::choose_pair(NeighRequest *rq)
{
  if (includegroup && rq->ghost)
    error->all(FLERR,"Neighbor include group not allowed with ghost neighbors");
  bool newtflag;
  if (rq->newton == 0 && newton_pair) newtflag = true;
  else if (rq->newton == 0 && !newton_pair) newtflag = false;
  else if (rq->newton == 1) newtflag = true;
  else if (rq->newton == 2) newtflag = false;
  else error->all(FLERR,"Illegal 'newton' flag in neighbor list request");
  int mask;
  for (int i = 0; i < npclass; i++) {
    mask = pairmasks[i];
    if (rq->copy) {
      if (!(mask & NP_COPY)) continue;
      if (rq->trim) {
        if (!rq->trim != !(mask & NP_TRIM)) continue;
        if (!rq->omp != !(mask & NP_OMP)) continue;
      }
      return i+1;
    }
    if (rq->half) {
      if (!(mask & NP_HALF)) continue;
    } else if (rq->full) {
      if (!(mask & NP_FULL)) continue;
    }
    if (newtflag) {
      if (!(mask & NP_NEWTON)) continue;
    } else if (!newtflag) {
      if (!(mask & NP_NEWTOFF)) continue;
    }
    if (mask & NP_MOLONLY) continue;
    if (!rq->ghost != !(mask & NP_GHOST)) continue;
    if (!rq->size != !(mask & NP_SIZE)) continue;
    if (!rq->respaouter != !(mask & NP_RESPA)) continue;
    if (!rq->granonesided != !(mask & NP_ONESIDE)) continue;
    if (!rq->respaouter != !(mask & NP_RESPA)) continue;
    if (!rq->bond != !(mask & NP_BOND)) continue;
    if (!rq->omp != !(mask & NP_OMP)) continue;
    if (!rq->ssa != !(mask & NP_SSA)) continue;
    if (!rq->skip != !(mask & NP_SKIP)) continue;
    if (!rq->trim != !(mask & NP_TRIM)) continue;
    if (!rq->halffull != !(mask & NP_HALF_FULL)) continue;
    if (!rq->off2on != !(mask & NP_OFF2ON)) continue;
    if (style == Neighbor::NSQ) {
      if (!(mask & NP_NSQ)) continue;
    } else if (style == Neighbor::BIN) {
      if (!(mask & NP_BIN)) continue;
    } else if (style == Neighbor::MULTI_OLD) {
      if (!(mask & NP_MULTI_OLD)) continue;
    } else if (style == Neighbor::MULTI) {
      if (!(mask & NP_MULTI)) continue;
    }
    if (triclinic) {
      if (!(mask & NP_TRI)) continue;
    } else if (!triclinic) {
      if (!(mask & NP_ORTHO)) continue;
    }
    return i+1;
  }
  return -1;
}
int Neighbor::request(void *requestor, int instance)
{
  if (nrequest == maxrequest) {
    maxrequest += RQDELTA;
    requests = (NeighRequest **)
      memory->srealloc(requests,maxrequest*sizeof(NeighRequest *), "neighbor:requests");
  }
  requests[nrequest] = new NeighRequest(lmp, requestor, instance);
  nrequest++;
  return nrequest-1;
}
NeighRequest *Neighbor::add_request(Pair *requestor, int flags)
{
  int irequest = request(requestor, requestor->instance_me);
  auto req = requests[irequest];
  req->apply_flags(flags);
  if (requestor->suffix_flag & Suffix::INTEL) {
    req->omp = 0;
  }
  return req;
}
NeighRequest *Neighbor::add_request(Fix *requestor, int flags)
{
  int irequest = request(requestor, requestor->instance_me);
  auto req = requests[irequest];
  req->pair = 0;
  req->fix = 1;
  req->apply_flags(flags);
  return req;
}
NeighRequest *Neighbor::add_request(Compute *requestor, int flags)
{
  int irequest = request(requestor, requestor->instance_me);
  auto req = requests[irequest];
  req->pair = 0;
  req->compute = 1;
  req->apply_flags(flags);
  return req;
}
NeighRequest *Neighbor::add_request(Command *requestor, const char *style, int flags)
{
  int irequest = request(requestor, 0);
  auto req = requests[irequest];
  req->pair = 0;
  req->command = 1;
  req->occasional = 1;
  req->command_style = style;
  req->apply_flags(flags);
  return req;
}
void Neighbor::setup_bins()
{
  for (int i = 0; i < nbin; i++)
    neigh_bin[i]->setup_bins(style);
  last_setup_bins = update->ntimestep;
}
int Neighbor::decide()
{
  ago++;
  if (ago >= delay && ago % every == 0) {
    if (build_once) return 0;
    if (dist_check == 0) return 1;
    return check_distance();
  } else return 0;
}
int Neighbor::check_distance()
{
  double delx,dely,delz,rsq;
  double delta,deltasq,delta1,delta2;
  if (boxcheck) {
    if (triclinic == 0) {
      delx = bboxlo[0] - boxlo_hold[0];
      dely = bboxlo[1] - boxlo_hold[1];
      delz = bboxlo[2] - boxlo_hold[2];
      delta1 = sqrt(delx*delx + dely*dely + delz*delz);
      delx = bboxhi[0] - boxhi_hold[0];
      dely = bboxhi[1] - boxhi_hold[1];
      delz = bboxhi[2] - boxhi_hold[2];
      delta2 = sqrt(delx*delx + dely*dely + delz*delz);
      delta = 0.5 * (skin - (delta1+delta2));
      if (delta < 0.0) delta = 0.0;
      deltasq = delta*delta;
    } else {
      domain->box_corners();
      delta1 = delta2 = 0.0;
      for (int i = 0; i < 8; i++) {
        delx = corners[i][0] - corners_hold[i][0];
        dely = corners[i][1] - corners_hold[i][1];
        delz = corners[i][2] - corners_hold[i][2];
        delta = sqrt(delx*delx + dely*dely + delz*delz);
        if (delta > delta1) delta1 = delta;
        else if (delta > delta2) delta2 = delta;
      }
      delta = 0.5 * (skin - (delta1+delta2));
      if (delta < 0.0) delta = 0.0;
      deltasq = delta*delta;
    }
  } else deltasq = triggersq;
  double **x = atom->x;
  int nlocal = atom->nlocal;
  if (includegroup) nlocal = atom->nfirst;
  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    delx = x[i][0] - xhold[i][0];
    dely = x[i][1] - xhold[i][1];
    delz = x[i][2] - xhold[i][2];
    rsq = delx*delx + dely*dely + delz*delz;
    if (rsq > deltasq) flag = 1;
  }
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_MAX,world);
  if (flagall && ago == MAX(every,delay)) ndanger++;
  return flagall;
}
void Neighbor::build(int topoflag)
{
  int i,m;
  ago = 0;
  ncalls++;
  lastcall = update->ntimestep;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  if (style == Neighbor::MULTI) build_collection(0);
  if (nall > NEIGHMASK)
    error->one(FLERR,"Too many local+ghost atoms for neighbor list");
  if (dist_check) {
    double **x = atom->x;
    if (includegroup) nlocal = atom->nfirst;
    if (atom->nmax > maxhold) {
      maxhold = atom->nmax;
      memory->destroy(xhold);
      memory->create(xhold,maxhold,3,"neigh:xhold");
    }
    for (i = 0; i < nlocal; i++) {
      xhold[i][0] = x[i][0];
      xhold[i][1] = x[i][1];
      xhold[i][2] = x[i][2];
    }
    if (boxcheck) {
      if (triclinic == 0) {
        boxlo_hold[0] = bboxlo[0];
        boxlo_hold[1] = bboxlo[1];
        boxlo_hold[2] = bboxlo[2];
        boxhi_hold[0] = bboxhi[0];
        boxhi_hold[1] = bboxhi[1];
        boxhi_hold[2] = bboxhi[2];
      } else {
        domain->box_corners();
        corners = domain->corners;
        for (i = 0; i < 8; i++) {
          corners_hold[i][0] = corners[i][0];
          corners_hold[i][1] = corners[i][1];
          corners_hold[i][2] = corners[i][2];
        }
      }
    }
  }
  if (style != Neighbor::NSQ) {
    if (last_setup_bins < 0) setup_bins();
    for (i = 0; i < nbin; i++) {
      neigh_bin[i]->bin_atoms_setup(nall);
      neigh_bin[i]->bin_atoms();
    }
  }
  for (i = 0; i < npair_perpetual; i++) {
    m = plist[i];
    if (!lists[m]->copy || lists[m]->trim)
      lists[m]->grow(nlocal,nall);
    neigh_pair[m]->build_setup();
    neigh_pair[m]->build(lists[m]);
  }
}
void Neighbor::build_one(class NeighList *mylist, int preflag)
{
  if (mylist == nullptr)
    error->all(FLERR,"Trying to build an occasional neighbor list before initialization complete");
  if (!mylist->occasional) error->all(FLERR,"Neighbor::build_one() invoked on perpetual list");
  NPair *np = neigh_pair[mylist->index];
  if (preflag) {
    if (np->last_build > lastcall) return;
  } else {
    if (np->last_build >= lastcall) return;
  }
  if (mylist->listcopy && mylist->listcopy->occasional)
    build_one(mylist->listcopy,preflag);
  if (mylist->listfull && mylist->listfull->occasional)
    build_one(mylist->listfull,preflag);
  if (mylist->listskip && mylist->listskip->occasional)
    build_one(mylist->listskip,preflag);
  if (!mylist->copy || mylist->trim)
    mylist->grow(atom->nlocal,atom->nlocal+atom->nghost);
  np->build_setup();
  np->build(mylist);
}
void Neighbor::set(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal neighbor command: expected 2 arguments but found {}", narg);
  skin = utils::numeric(FLERR,arg[0],false,lmp);
  if (skin < 0.0) error->all(FLERR, "Invalid neighbor argument: {}", arg[0]);
  if (strcmp(arg[1],"nsq") == 0) style = Neighbor::NSQ;
  else if (strcmp(arg[1],"bin") == 0) style = Neighbor::BIN;
  else if (strcmp(arg[1],"multi") == 0) {
    style = Neighbor::MULTI;
    ncollections = atom->ntypes;
  } else if (strcmp(arg[1],"multi/old") == 0) style = Neighbor::MULTI_OLD;
  else error->all(FLERR,"Unknown neighbor {} argument: {}", arg[0], arg[1]);
}
void Neighbor::reset_timestep(bigint )
{
  for (int i = 0; i < nbin; i++)
    neigh_bin[i]->last_bin = -1;
  for (int i = 0; i < nlist; i++) {
    if (!neigh_pair[i]) continue;
    neigh_pair[i]->last_build = -1;
  }
  lastcall = -1;
  last_setup_bins = -1;
}
void Neighbor::modify_params(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "neigh_modify every", error);
      every = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (every <= 0) error->all(FLERR, "Invalid neigh_modify every argument: {}", every);
      iarg += 2;
    } else if (strcmp(arg[iarg],"delay") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "neigh_modify delay", error);
      delay = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (delay < 0) error->all(FLERR, "Invalid neigh_modify delay argument: {}", delay);
      iarg += 2;
    } else if (strcmp(arg[iarg],"check") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "neigh_modify check", error);
      dist_check = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"once") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "neigh_modify once", error);
      build_once = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"page") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "neigh_modify page", error);
      old_pgsize = pgsize;
      pgsize = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"one") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "neigh_modify one", error);
      old_oneatom = oneatom;
      oneatom = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"binsize") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "neigh_modify binsize", error);
      binsize_user = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (binsize_user <= 0.0) binsizeflag = 0;
      else binsizeflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"cluster") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "neigh_modify cluster", error);
      cluster_check = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"include") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "neigh_modify include", error);
      includegroup = group->find(arg[iarg+1]);
      if (includegroup < 0)
        error->all(FLERR, "Invalid include keyword: group {} not found", arg[iarg+1]);
      if (atom->firstgroupname == nullptr)
          error->all(FLERR, "Invalid include keyword: atom_modify first command must be used");
      if (strcmp(arg[iarg+1],atom->firstgroupname) != 0)
        error->all(FLERR, "Neigh_modify include group != atom_modify first group: {}", atom->firstgroupname);
      iarg += 2;
    } else if (strcmp(arg[iarg],"exclude") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "neigh_modify exclude", error);
      if (strcmp(arg[iarg+1],"type") == 0) {
        if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "neigh_modify exclude type", error);
        if (nex_type == maxex_type) {
          maxex_type += EXDELTA;
          memory->grow(ex1_type,maxex_type,"neigh:ex1_type");
          memory->grow(ex2_type,maxex_type,"neigh:ex2_type");
        }
        ex1_type[nex_type] = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
        ex2_type[nex_type] = utils::inumeric(FLERR,arg[iarg+3],false,lmp);
        nex_type++;
        iarg += 4;
      } else if (strcmp(arg[iarg+1],"group") == 0) {
        if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "neigh_modify exclude group", error);
        if (nex_group == maxex_group) {
          maxex_group += EXDELTA;
          memory->grow(ex1_group,maxex_group,"neigh:ex1_group");
          memory->grow(ex2_group,maxex_group,"neigh:ex2_group");
        }
        ex1_group[nex_group] = group->find(arg[iarg+2]);
        ex2_group[nex_group] = group->find(arg[iarg+3]);
        if (ex1_group[nex_group] == -1)
          error->all(FLERR, "Invalid exclude group keyword: group {} not found", arg[iarg+2]);
        if (ex2_group[nex_group] == -1)
            error->all(FLERR, "Invalid exclude group keyword: group {} not found", arg[iarg+3]);
        nex_group++;
        iarg += 4;
      } else if (strcmp(arg[iarg+1],"molecule/inter") == 0 ||
                 strcmp(arg[iarg+1],"molecule/intra") == 0) {
        if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "neigh_modify exclude molecule", error);
        if (nex_mol == maxex_mol) {
          maxex_mol += EXDELTA;
          memory->grow(ex_mol_group,maxex_mol,"neigh:ex_mol_group");
   memory->grow(ex_mol_intra,maxex_mol,"neigh:ex_mol_intra");
        }
        ex_mol_group[nex_mol] = group->find(arg[iarg+2]);
        if (ex_mol_group[nex_mol] == -1)
          error->all(FLERR, "Invalid exclude keyword:group {} not found", arg[iarg+2]);
        if (strcmp(arg[iarg+1],"molecule/intra") == 0)
          ex_mol_intra[nex_mol] = 1;
        else
          ex_mol_intra[nex_mol] = 0;
        nex_mol++;
        iarg += 3;
      } else if (strcmp(arg[iarg+1],"none") == 0) {
        nex_type = nex_group = nex_mol = 0;
        iarg += 2;
      } else error->all(FLERR,"Unknown neigh_modify exclude keyword: {}", arg[iarg+1]);
    } else if (strcmp(arg[iarg],"collection/interval") == 0) {
      if (style != Neighbor::MULTI)
        error->all(FLERR,"Cannot use collection/interval command without multi setting");
      if (iarg+2 > narg)
        utils::missing_cmd_args(FLERR, "neigh_modify collection/interval", error);
      ncollections = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (ncollections < 1)
        error->all(FLERR, "Invalid collection/interval keyword: illegal number of custom collections: {}", ncollections);
      if (iarg+2+ncollections > narg)
        error->all(FLERR, "Invalid collection/interval keyword: expected {} separate lists of types", ncollections);
      int i;
      comm->ncollections_cutoff = 0;
      interval_collection_flag = 1;
      custom_collection_flag = 1;
      memory->grow(collection2cut,ncollections,"neigh:collection2cut");
      double cut_interval;
      for (i = 0; i < ncollections; i++){
        cut_interval = utils::numeric(FLERR,arg[iarg+2+i],false,lmp);
        collection2cut[i] = cut_interval;
        if (i != 0)
          if (collection2cut[i-1] >= collection2cut[i])
            error->all(FLERR,"Nonsequential interval cutoffs in collection/interval setting");
      }
      iarg += 2 + ncollections;
    } else if (strcmp(arg[iarg],"collection/type") == 0) {
      if (style != Neighbor::MULTI)
        error->all(FLERR,"Cannot use collection/type command without multi setting");
      if (iarg+2 > narg)
        utils::missing_cmd_args(FLERR, "neigh_modify collection/type", error);
      ncollections = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (ncollections < 1)
        error->all(FLERR, "Invalid collection/type keyword: illegal number of custom collections: {}", ncollections);
      if (iarg+2+ncollections > narg)
        error->all(FLERR, "Invalid collection/type keyword: expected {} separate lists of types", ncollections);
      int ntypes = atom->ntypes;
      int nlo, nhi, i, k;
      comm->ncollections_cutoff = 0;
      interval_collection_flag = 0;
      custom_collection_flag = 1;
      if (!type2collection)
        memory->create(type2collection,ntypes+1,"neigh:type2collection");
      for (i = 1; i <= ntypes; i++)
        type2collection[i] = -1;
      for (i = 0; i < ncollections; i++){
        std::vector<std::string> words = Tokenizer(arg[iarg+2+i], ",").as_vector();
        for (const auto &word : words) {
          utils::bounds(FLERR,word,1,ntypes,nlo,nhi,error);
          for (k = nlo; k <= nhi; k++) {
            if (type2collection[k] != -1)
              error->all(FLERR,"Type specified more than once in collection/type commnd");
            type2collection[k] = i;
          }
        }
      }
      for (i = 1; i <= ntypes; i++){
        if (type2collection[i] == -1) {
          error->all(FLERR,"Type missing in collection/type commnd");
        }
      }
      iarg += 2 + ncollections;
    } else error->all(FLERR,"Unknown neigh_modify keyword: {}", arg[iarg]);
  }
}
void Neighbor::modify_params(const std::string &modcmd)
{
  auto args = utils::split_words(modcmd);
  auto newarg = new char*[args.size()];
  int i=0;
  for (const auto &arg : args) {
    newarg[i++] = (char *)arg.c_str();
  }
  modify_params(args.size(),newarg);
  delete[] newarg;
}
void Neighbor::exclusion_group_group_delete(int group1, int group2)
{
  int m, mlast;
  for (m = 0; m < nex_group; m++)
    if (ex1_group[m] == group1 && ex2_group[m] == group2 )
      break;
  mlast = m;
  if (mlast == nex_group)
    error->all(FLERR,"Unable to find group-group exclusion");
  for (m = mlast+1; m < nex_group; m++) {
    ex1_group[m-1] = ex1_group[m];
    ex2_group[m-1] = ex2_group[m];
    ex1_bit[m-1] = ex1_bit[m];
    ex2_bit[m-1] = ex2_bit[m];
  }
  nex_group--;
}
int Neighbor::exclude_setting()
{
  return exclude;
}
void Neighbor::set_overlap_topo(int s)
{
  overlap_topo = s;
}
int Neighbor::any_full()
{
  int any_full = 0;
  for (int i = 0; i < old_nrequest; i++) {
    if (old_requests[i]->full) any_full = 1;
  }
  return any_full;
}
void Neighbor::build_collection(int istart)
{
  if (style != Neighbor::MULTI)
    error->all(FLERR, "Cannot define atom collections without neighbor style multi");
  int nmax = atom->nlocal+atom->nghost;
  if (nmax > nmax_collection) {
    nmax_collection = nmax+DELTA_PERATOM;
    memory->grow(collection, nmax_collection, "neigh:collection");
  }
  if (finite_cut_flag) {
    double cut;
    int icollection;
    for (int i = istart; i < nmax; i++){
      cut = force->pair->atom2cut(i);
      collection[i] = -1;
      for (icollection = 0; icollection < ncollections; icollection++){
        if (collection2cut[icollection] >= cut) {
          collection[i] = icollection;
          break;
        }
      }
      if (collection[i] == -1)
        error->one(FLERR, "Atom cutoff exceeds interval cutoffs for multi");
    }
  } else {
    int *type = atom->type;
    for (int i = istart; i < nmax; i++){
      collection[i] = type2collection[type[i]];
    }
  }
}
bigint Neighbor::get_nneigh_full()
{
  int m;
  for (m = 0; m < old_nrequest; m++)
    if (old_requests[m]->full && !old_requests[m]->skip) break;
  bigint nneighfull = -1;
  if (m < old_nrequest) {
    nneighfull = 0;
    if (lists[m]->numneigh) {
      int inum = neighbor->lists[m]->inum;
      int *ilist = neighbor->lists[m]->ilist;
      int *numneigh = neighbor->lists[m]->numneigh;
      for (int i = 0; i < inum; i++)
        nneighfull += numneigh[ilist[i]];
    };
  }
  return nneighfull;
}
bigint Neighbor::get_nneigh_half()
{
  int m;
  for (m = 0; m < old_nrequest; m++)
    if (old_requests[m]->half && !old_requests[m]->skip && lists[m] && lists[m]->numneigh) break;
  bigint nneighhalf = -1;
  if (m < old_nrequest) {
    nneighhalf = 0;
  }
  return nneighhalf;
}
double Neighbor::memory_usage()
{
  double bytes = 0;
  bytes += memory->usage(xhold,maxhold,3);
  for (int i = 0; i < nlist; i++)
    if (lists[i]) bytes += lists[i]->memory_usage();
  for (int i = 0; i < nbin; i++)
    bytes += neigh_bin[i]->memory_usage();
  return bytes;
}
