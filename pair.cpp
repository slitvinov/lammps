#include "pair.h"
#include "atom.h"
#include "atom_masks.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neighbor.h"
#include "update.h"
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstring>
using namespace LAMMPS_NS;
using MathConst::MY_ISPI4;
using MathConst::THIRD;
enum { NONE, RLINEAR, RSQ, BMP };
static const std::string mixing_rule_names[Pair::SIXTHPOWER + 1] = {
    "geometric", "arithmetic", "sixthpower"};
int Pair::instance_total = 0;
Pair::Pair(LAMMPS *lmp)
    : Pointers(lmp), eatom(nullptr), vatom(nullptr), cvatom(nullptr),
      cutsq(nullptr), setflag(nullptr), cutghost(nullptr), rtable(nullptr),
      drtable(nullptr), ftable(nullptr), dftable(nullptr), ctable(nullptr),
      dctable(nullptr), etable(nullptr), detable(nullptr), ptable(nullptr),
      dptable(nullptr), vtable(nullptr), dvtable(nullptr), rdisptable(nullptr),
      drdisptable(nullptr), fdisptable(nullptr), dfdisptable(nullptr),
      edisptable(nullptr), dedisptable(nullptr), pvector(nullptr),
      svector(nullptr), list(nullptr), listhalf(nullptr), listfull(nullptr),
      list_tally_compute(nullptr), elements(nullptr), elem1param(nullptr),
      elem2param(nullptr), elem3param(nullptr), map(nullptr) {
  instance_me = instance_total++;
  eng_vdwl = eng_coul = 0.0;
  comm_forward = comm_reverse = comm_reverse_off = 0;
  single_enable = 1;
  born_matrix_enable = 0;
  restartinfo = 1;
  respa_enable = 0;
  one_coeff = 0;
  no_virial_fdotr_compute = 0;
  writedata = 0;
  finitecutflag = 0;
  ghostneigh = 0;
  unit_convert_flag = utils::NOCONVERT;
  did_mix = false;
  nextra = 0;
  single_extra = 0;
  ewaldflag = pppmflag = msmflag = dispersionflag = tip4pflag = dipoleflag =
      spinflag = 0;
  reinitflag = 1;
  centroidstressflag = CENTROID_SAME;
  compute_flag = 1;
  manybody_flag = 0;
  offset_flag = 0;
  mix_flag = GEOMETRIC;
  mixed_flag = 1;
  tail_flag = 0;
  etail = ptail = etail_ij = ptail_ij = 0.0;
  ncoultablebits = 12;
  ndisptablebits = 12;
  tabinner = sqrt(2.0);
  tabinner_disp = sqrt(2.0);
  trim_flag = 1;
  allocated = 0;
  maxeatom = maxvatom = maxcvatom = 0;
  num_tally_compute = 0;
  nelements = nparams = maxparam = 0;
  nondefault_history_transfer = 0;
  beyond_contact = 0;
  execution_space = Host;
  datamask_read = ALL_MASK;
  datamask_modify = ALL_MASK;
  reverse_comm_device = 0;
  copymode = 0;
}
Pair::~Pair() {
  num_tally_compute = 0;
  memory->sfree((void *)list_tally_compute);
  list_tally_compute = nullptr;
  if (copymode)
    return;
  if (elements)
    for (int i = 0; i < nelements; i++)
      delete[] elements[i];
  delete[] elements;
  delete[] map;
  memory->destroy(eatom);
  memory->destroy(vatom);
  memory->destroy(cvatom);
}
void Pair::init() {
  int i, j;
  if (offset_flag && tail_flag)
    error->all(FLERR, "Cannot have both pair_modify shift and tail set to yes");
  if (tail_flag && domain->dimension == 2)
    error->all(FLERR, "Cannot use pair tail corrections with 2d simulations");
  if (tail_flag && domain->nonperiodic && comm->me == 0)
    error->warning(FLERR,
                   "Using pair tail corrections with non-periodic system");
  if (!compute_flag && tail_flag && comm->me == 0)
    error->warning(FLERR,
                   "Using pair tail corrections with pair_modify compute no");
  if (!compute_flag && offset_flag && comm->me == 0)
    error->warning(FLERR,
                   "Using pair potential shift with pair_modify compute no");
  if (!allocated)
    error->all(FLERR, "All pair coeffs are not set");
  for (i = 1; i <= atom->ntypes; i++)
    if (setflag[i][i] == 0)
      error->all(FLERR, "All pair coeffs are not set");
  init_style();
  cutforce = 0.0;
  etail = ptail = 0.0;
  mixed_flag = 1;
  double cut;
  int mixed_count = 0;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      did_mix = false;
      cut = init_one(i, j);
      cutsq[i][j] = cutsq[j][i] = cut * cut;
      cutforce = MAX(cutforce, cut);
      if (i != j) {
        if (setflag[i][j])
          mixed_flag = 0;
        if (did_mix)
          ++mixed_count;
      }
      if (tail_flag) {
        etail += etail_ij;
        ptail += ptail_ij;
        if (i != j) {
          etail += etail_ij;
          ptail += ptail_ij;
        }
      }
    }
  if (!manybody_flag && (comm->me == 0)) {
    const int num_mixed_pairs = atom->ntypes * (atom->ntypes - 1) / 2;
    if (utils::strmatch(force->pair_style, "^lj/class2"))
      utils::logmesg(
          lmp,
          "Generated {} of {} mixed pair_coeff terms from {}/{} mixing rule\n",
          mixed_count, num_mixed_pairs, "sixthpower",
          mixing_rule_names[mix_flag]);
    else
      utils::logmesg(
          lmp,
          "Generated {} of {} mixed pair_coeff terms from {} mixing rule\n",
          mixed_count, num_mixed_pairs, mixing_rule_names[mix_flag]);
  }
}
void Pair::init_style() { neighbor->add_request(this); }
void Pair::init_list(int, NeighList *ptr) { list = ptr; }
void Pair::compute_dummy(int eflag, int vflag) { ev_init(eflag, vflag); }
void Pair::add_tally_callback(Compute *ptr) {
  int i, found = -1;
  for (i = 0; i < num_tally_compute; ++i) {
    if (list_tally_compute[i] == ptr)
      found = i;
  }
  if (found < 0) {
    found = num_tally_compute;
    ++num_tally_compute;
    void *p = memory->srealloc((void *)list_tally_compute,
                               sizeof(Compute *) * num_tally_compute,
                               "pair:list_tally_compute");
    list_tally_compute = (Compute **)p;
    list_tally_compute[num_tally_compute - 1] = ptr;
  }
}
void Pair::del_tally_callback(Compute *ptr) {
  int i, found = -1;
  for (i = 0; i < num_tally_compute; ++i) {
    if (list_tally_compute[i] == ptr)
      found = i;
  }
  if (found < 0)
    return;
  --num_tally_compute;
  for (i = found; i < num_tally_compute; ++i) {
    list_tally_compute[i] = list_tally_compute[i + 1];
  }
}
void Pair::map_element2type(int narg, char **arg, bool update_setflag) {
  int i, j;
  const int ntypes = atom->ntypes;
  if (narg != ntypes)
    error->all(FLERR, "Number of element to type mappings does not match "
                      "number of atom types");
  if (elements) {
    for (i = 0; i < nelements; i++)
      delete[] elements[i];
    delete[] elements;
  }
  elements = new char *[ntypes];
  for (i = 0; i < ntypes; i++)
    elements[i] = nullptr;
  nelements = 0;
  map[0] = -1;
  for (i = 1; i <= narg; i++) {
    std::string entry = arg[i - 1];
    if (entry == "NULL") {
      map[i] = -1;
      continue;
    }
    for (j = 0; j < nelements; j++)
      if (entry == elements[j])
        break;
    map[i] = j;
    if (j == nelements) {
      elements[j] = utils::strdup(entry);
      nelements++;
    }
  }
  if (update_setflag) {
    int count = 0;
    for (i = 1; i <= ntypes; i++) {
      for (j = i; j <= ntypes; j++) {
        setflag[i][j] = 0;
        if ((map[i] >= 0) && (map[j] >= 0)) {
          setflag[i][j] = 1;
          count++;
        }
      }
    }
    if (count == 0)
      error->all(FLERR, "Incorrect args for pair coefficients");
  }
}
void Pair::ev_setup(int eflag, int vflag, int alloc) {
  int i, n;
  eflag_either = eflag;
  eflag_global = eflag & ENERGY_GLOBAL;
  eflag_atom = eflag & ENERGY_ATOM;
  vflag_global = vflag & VIRIAL_PAIR;
  if (vflag & VIRIAL_FDOTR && no_virial_fdotr_compute == 1)
    vflag_global = 1;
  vflag_fdotr = 0;
  if (vflag & VIRIAL_FDOTR && no_virial_fdotr_compute == 0)
    vflag_fdotr = 1;
  vflag_atom = vflag & VIRIAL_ATOM;
  if (vflag & VIRIAL_CENTROID && centroidstressflag != CENTROID_AVAIL)
    vflag_atom = 1;
  cvflag_atom = 0;
  if (vflag & VIRIAL_CENTROID && centroidstressflag == CENTROID_AVAIL)
    cvflag_atom = 1;
  vflag_either = vflag_global || vflag_atom || cvflag_atom;
  evflag = eflag_either || vflag_either;
  if (eflag_atom && atom->nmax > maxeatom) {
    maxeatom = atom->nmax;
    if (alloc) {
      memory->destroy(eatom);
      memory->create(eatom, comm->nthreads * maxeatom, "pair:eatom");
    }
  }
  if (vflag_atom && atom->nmax > maxvatom) {
    maxvatom = atom->nmax;
    if (alloc) {
      memory->destroy(vatom);
      memory->create(vatom, comm->nthreads * maxvatom, 6, "pair:vatom");
    }
  }
  if (cvflag_atom && atom->nmax > maxcvatom) {
    maxcvatom = atom->nmax;
    if (alloc) {
      memory->destroy(cvatom);
      memory->create(cvatom, comm->nthreads * maxcvatom, 9, "pair:cvatom");
    }
  }
  if (eflag_global)
    eng_vdwl = eng_coul = 0.0;
  if (vflag_global || vflag_fdotr)
    for (i = 0; i < 6; i++)
      virial[i] = 0.0;
  if (eflag_atom && alloc) {
    n = atom->nlocal;
    if (force->newton)
      n += atom->nghost;
    for (i = 0; i < n; i++)
      eatom[i] = 0.0;
  }
  if (vflag_atom && alloc) {
    n = atom->nlocal;
    if (force->newton)
      n += atom->nghost;
    for (i = 0; i < n; i++) {
      vatom[i][0] = 0.0;
      vatom[i][1] = 0.0;
      vatom[i][2] = 0.0;
      vatom[i][3] = 0.0;
      vatom[i][4] = 0.0;
      vatom[i][5] = 0.0;
    }
  }
  if (cvflag_atom && alloc) {
    n = atom->nlocal;
    if (force->newton)
      n += atom->nghost;
    for (i = 0; i < n; i++) {
      cvatom[i][0] = 0.0;
      cvatom[i][1] = 0.0;
      cvatom[i][2] = 0.0;
      cvatom[i][3] = 0.0;
      cvatom[i][4] = 0.0;
      cvatom[i][5] = 0.0;
      cvatom[i][6] = 0.0;
      cvatom[i][7] = 0.0;
      cvatom[i][8] = 0.0;
      cvatom[i][9] = 0.0;
    }
  }
  if (num_tally_compute > 0) {
    for (int k = 0; k < num_tally_compute; ++k) {
      Compute *c = list_tally_compute[k];
      c->pair_setup_callback(eflag, vflag);
    }
  }
}
void Pair::ev_unset() {
  evflag = 0;
  eflag_either = 0;
  eflag_global = 0;
  eflag_atom = 0;
  vflag_either = 0;
  vflag_global = 0;
  vflag_atom = 0;
  cvflag_atom = 0;
  vflag_fdotr = 0;
}
