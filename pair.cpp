#include <map>
#include <set>
#include <unordered_set>
#include <vector>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstring>
#include <string>
#include "utils.h"
#include "lammps.h"
#include "pointers.h"
#include "pair.h"
#include "atom.h"
#include "atom_masks.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neighbor.h"
#include "update.h"
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
void Pair::init() {
  int i, j;
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
    }
}
void Pair::init_style() { neighbor->add_request(this); }
void Pair::init_list(int, NeighList *ptr) { list = ptr; }
