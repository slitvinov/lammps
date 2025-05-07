#include <set>
#include <map>
#include <vector>
#include "utils.h"
#include "fix.h"
#include "atom.h"
#include "atom_masks.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include <cstring>
using namespace LAMMPS_NS;
using namespace FixConst;
int Fix::instance_total = 0;
Fix::Fix(LAMMPS *lmp, int, char **arg)
    : Pointers(lmp), id(nullptr), style(nullptr), extlist(nullptr),
      vector_atom(nullptr), array_atom(nullptr), vector_local(nullptr),
      array_local(nullptr), eatom(nullptr), vatom(nullptr), cvatom(nullptr) {
  instance_me = instance_total++;
  id = utils::strdup(arg[0]);
  igroup = group->find(arg[1]);
  groupbit = group->bitmask[igroup];
  style = utils::strdup(arg[2]);
  restart_global = restart_peratom = restart_file = 0;
  force_reneighbor = 0;
  box_change = NO_BOX_CHANGE;
  thermo_energy = 0;
  thermo_virial = 0;
  energy_global_flag = energy_peratom_flag = 0;
  virial_global_flag = virial_peratom_flag = 0;
  ecouple_flag = 0;
  rigid_flag = 0;
  no_change_box = 0;
  time_integrate = 0;
  time_depend = 0;
  create_attribute = 0;
  restart_pbc = 0;
  wd_header = wd_section = 0;
  dynamic_group_allow = 0;
  dynamic = 0;
  dof_flag = 0;
  special_alter_flag = 0;
  enforce2d_flag = 0;
  respa_level_support = 0;
  respa_level = -1;
  maxexchange = 0;
  maxexchange_dynamic = 0;
  pre_exchange_migrate = 0;
  stores_ids = 0;
  scalar_flag = vector_flag = array_flag = 0;
  peratom_flag = local_flag = pergrid_flag = 0;
  global_freq = local_freq = peratom_freq = pergrid_freq = -1;
  size_vector_variable = size_array_rows_variable = 0;
  comm_forward = comm_reverse = comm_border = 0;
  restart_reset = 0;
  nevery = 1;
  global_freq = 1;
  maxeatom = maxvatom = maxcvatom = 0;
  vflag_atom = cvflag_atom = 0;
  centroidstressflag = CENTROID_SAME;
  execution_space = Host;
  datamask_read = ALL_MASK;
  datamask_modify = ALL_MASK;
  forward_comm_device = 0;
  copymode = 0;
}
Fix::~Fix() {
  if (copymode)
    return;
  delete[] id;
  delete[] style;
  memory->destroy(eatom);
  memory->destroy(vatom);
  memory->destroy(cvatom);
}
void Fix::modify_params(int narg, char **arg) {
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "dynamic/dof") == 0) {
      dynamic = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "energy") == 0) {
      thermo_energy = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "virial") == 0) {
      thermo_virial = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "respa") == 0) {
      int lvl = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      respa_level = lvl - 1;
      iarg += 2;
    } else {
      int n = modify_param(narg - iarg, &arg[iarg]);
      iarg += n;
    }
  }
}
void Fix::ev_setup(int eflag, int vflag) {
  int i, n;
  evflag = 1;
  if (!thermo_energy)
    eflag_either = eflag_global = eflag_atom = 0;
  else {
    eflag_either = eflag;
    eflag_global = eflag & ENERGY_GLOBAL;
    eflag_atom = eflag & ENERGY_ATOM;
  }
  if (!thermo_virial)
    vflag_either = vflag_global = vflag_atom = 0;
  else {
    vflag_either = vflag;
    vflag_global = vflag & (VIRIAL_PAIR | VIRIAL_FDOTR);
    if (centroidstressflag != CENTROID_AVAIL) {
      vflag_atom = vflag & (VIRIAL_ATOM | VIRIAL_CENTROID);
      cvflag_atom = 0;
    } else {
      vflag_atom = vflag & VIRIAL_ATOM;
      cvflag_atom = vflag & VIRIAL_CENTROID;
    }
  }
  if (eflag_atom && atom->nlocal > maxeatom) {
    maxeatom = atom->nmax;
    memory->destroy(eatom);
    memory->create(eatom, maxeatom, "fix:eatom");
  }
  if (vflag_atom && atom->nlocal > maxvatom) {
    maxvatom = atom->nmax;
    memory->destroy(vatom);
    memory->create(vatom, maxvatom, 6, "fix:vatom");
  }
  if (cvflag_atom && atom->nlocal > maxcvatom) {
    maxcvatom = atom->nmax;
    memory->destroy(cvatom);
    memory->create(cvatom, maxcvatom, 9, "fix:cvatom");
  }
  if (vflag_global)
    for (i = 0; i < 6; i++)
      virial[i] = 0.0;
  if (eflag_atom) {
    n = atom->nlocal;
    for (i = 0; i < n; i++)
      eatom[i] = 0.0;
  }
  if (vflag_atom) {
    n = atom->nlocal;
    for (i = 0; i < n; i++) {
      vatom[i][0] = 0.0;
      vatom[i][1] = 0.0;
      vatom[i][2] = 0.0;
      vatom[i][3] = 0.0;
      vatom[i][4] = 0.0;
      vatom[i][5] = 0.0;
    }
  }
  if (cvflag_atom) {
    n = atom->nlocal;
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
    }
  }
}
void Fix::v_setup(int vflag) {
  int i, n;
  evflag = 1;
  vflag_global = vflag & (VIRIAL_PAIR | VIRIAL_FDOTR);
  if (centroidstressflag != CENTROID_AVAIL) {
    vflag_atom = vflag & (VIRIAL_ATOM | VIRIAL_CENTROID);
    cvflag_atom = 0;
  } else {
    vflag_atom = vflag & VIRIAL_ATOM;
    cvflag_atom = vflag & VIRIAL_CENTROID;
  }
  if (vflag_atom && atom->nlocal > maxvatom) {
    maxvatom = atom->nmax;
    memory->destroy(vatom);
    memory->create(vatom, maxvatom, 6, "fix:vatom");
  }
  if (cvflag_atom && atom->nlocal > maxcvatom) {
    maxcvatom = atom->nmax;
    memory->destroy(cvatom);
    memory->create(cvatom, maxcvatom, 9, "fix:cvatom");
  }
  if (vflag_global)
    for (i = 0; i < 6; i++)
      virial[i] = 0.0;
  if (vflag_atom) {
    n = atom->nlocal;
    for (i = 0; i < n; i++) {
      vatom[i][0] = 0.0;
      vatom[i][1] = 0.0;
      vatom[i][2] = 0.0;
      vatom[i][3] = 0.0;
      vatom[i][4] = 0.0;
      vatom[i][5] = 0.0;
    }
  }
  if (cvflag_atom) {
    n = atom->nlocal;
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
    }
  }
}
void Fix::ev_tally(int n, int *list, double total, double eng, double *v) {
  if (eflag_atom) {
    double fraction = eng / total;
    for (int i = 0; i < n; i++)
      eatom[list[i]] += fraction;
  }
  v_tally(n, list, total, v);
}
void Fix::v_tally(int n, int *list, double total, double *v) {
  int m;
  if (vflag_global) {
    double fraction = n / total;
    virial[0] += fraction * v[0];
    virial[1] += fraction * v[1];
    virial[2] += fraction * v[2];
    virial[3] += fraction * v[3];
    virial[4] += fraction * v[4];
    virial[5] += fraction * v[5];
  }
  if (vflag_atom) {
    double fraction = 1.0 / total;
    for (int i = 0; i < n; i++) {
      m = list[i];
      vatom[m][0] += fraction * v[0];
      vatom[m][1] += fraction * v[1];
      vatom[m][2] += fraction * v[2];
      vatom[m][3] += fraction * v[3];
      vatom[m][4] += fraction * v[4];
      vatom[m][5] += fraction * v[5];
    }
  }
}
void Fix::v_tally(int n, int *list, double total, double *vtot,
                  double rlist[][3], double flist[][3], double center[]) {
  v_tally(n, list, total, vtot);
  if (cvflag_atom) {
    for (int i = 0; i < n; i++) {
      const double ri0[3] = {
          rlist[i][0] - center[0],
          rlist[i][1] - center[1],
          rlist[i][2] - center[2],
      };
      cvatom[list[i]][0] += ri0[0] * flist[i][0];
      cvatom[list[i]][1] += ri0[1] * flist[i][1];
      cvatom[list[i]][2] += ri0[2] * flist[i][2];
      cvatom[list[i]][3] += ri0[0] * flist[i][1];
      cvatom[list[i]][4] += ri0[0] * flist[i][2];
      cvatom[list[i]][5] += ri0[1] * flist[i][2];
      cvatom[list[i]][6] += ri0[1] * flist[i][0];
      cvatom[list[i]][7] += ri0[2] * flist[i][0];
      cvatom[list[i]][8] += ri0[2] * flist[i][1];
    }
  }
}
void Fix::v_tally(int n, int *list, double total, double *vtot, int nlocal,
                  int npair, int pairlist[][2], double *fpairlist,
                  double dellist[][3]) {
  v_tally(n, list, total, vtot);
  if (cvflag_atom) {
    double v[6];
    for (int i = 0; i < npair; i++) {
      v[0] = 0.5 * dellist[i][0] * dellist[i][0] * fpairlist[i];
      v[1] = 0.5 * dellist[i][1] * dellist[i][1] * fpairlist[i];
      v[2] = 0.5 * dellist[i][2] * dellist[i][2] * fpairlist[i];
      v[3] = 0.5 * dellist[i][0] * dellist[i][1] * fpairlist[i];
      v[4] = 0.5 * dellist[i][0] * dellist[i][2] * fpairlist[i];
      v[5] = 0.5 * dellist[i][1] * dellist[i][2] * fpairlist[i];
      const int i0 = pairlist[i][0];
      const int i1 = pairlist[i][1];
      if (i0 < nlocal) {
        cvatom[i0][0] += v[0];
        cvatom[i0][1] += v[1];
        cvatom[i0][2] += v[2];
        cvatom[i0][3] += v[3];
        cvatom[i0][4] += v[4];
        cvatom[i0][5] += v[5];
        cvatom[i0][6] += v[3];
        cvatom[i0][7] += v[4];
        cvatom[i0][8] += v[5];
      }
      if (i1 < nlocal) {
        cvatom[i1][0] += v[0];
        cvatom[i1][1] += v[1];
        cvatom[i1][2] += v[2];
        cvatom[i1][3] += v[3];
        cvatom[i1][4] += v[4];
        cvatom[i1][5] += v[5];
        cvatom[i1][6] += v[3];
        cvatom[i1][7] += v[4];
        cvatom[i1][8] += v[5];
      }
    }
  }
}
void Fix::v_tally(int i, double *v) {
  if (vflag_global) {
    virial[0] += v[0];
    virial[1] += v[1];
    virial[2] += v[2];
    virial[3] += v[3];
    virial[4] += v[4];
    virial[5] += v[5];
  }
  if (vflag_atom) {
    vatom[i][0] += v[0];
    vatom[i][1] += v[1];
    vatom[i][2] += v[2];
    vatom[i][3] += v[3];
    vatom[i][4] += v[4];
    vatom[i][5] += v[5];
  }
}
void Fix::v_tally(int n, int i, double vn) {
  if (vflag_global)
    virial[n] += vn;
  if (vflag_atom)
    vatom[i][n] += vn;
}
