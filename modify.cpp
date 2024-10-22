#include "modify.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "group.h"
#include "input.h"
#include "memory.h"
#include "region.h"
#include "style_fix.h"
#include "update.h"
#include <cstring>
using namespace LAMMPS_NS;
using namespace FixConst;
#define DELTA 4
#define BIG 1.0e20
template <typename S, typename T>
static S *style_creator(LAMMPS *lmp, int narg, char **arg) {
  return new T(lmp, narg, arg);
}
Modify::Modify(LAMMPS *lmp) : Pointers(lmp) {
  nfix = maxfix = 0;
  n_initial_integrate = n_post_integrate = 0;
  n_pre_exchange = n_pre_neighbor = n_post_neighbor = 0;
  n_pre_force = n_pre_reverse = n_post_force_any = 0;
  n_final_integrate = n_end_of_step = 0;
  n_energy_couple = n_energy_global = n_energy_atom = 0;
  n_initial_integrate_respa = n_post_integrate_respa = 0;
  n_pre_force_respa = n_post_force_respa_any = n_final_integrate_respa = 0;
  n_min_pre_exchange = n_min_pre_force = n_min_pre_reverse = 0;
  n_min_post_force = n_min_energy = 0;
  n_timeflag = -1;
  fix = nullptr;
  fmask = nullptr;
  list_initial_integrate = list_post_integrate = nullptr;
  list_pre_exchange = list_pre_neighbor = list_post_neighbor = nullptr;
  list_pre_force = list_pre_reverse = nullptr;
  list_post_force = list_post_force_group = nullptr;
  list_final_integrate = list_end_of_step = nullptr;
  list_energy_couple = list_energy_global = list_energy_atom = nullptr;
  list_initial_integrate_respa = list_post_integrate_respa = nullptr;
  list_pre_force_respa = list_post_force_respa = nullptr;
  list_final_integrate_respa = nullptr;
  list_min_pre_exchange = list_min_pre_neighbor = list_min_post_neighbor =
      nullptr;
  list_min_pre_force = list_min_pre_reverse = list_min_post_force = nullptr;
  list_min_energy = nullptr;
  end_of_step_every = nullptr;
  list_timeflag = nullptr;
  restart_pbc_any = 0;
  nfix_restart_global = 0;
  id_restart_global = style_restart_global = nullptr;
  state_restart_global = nullptr;
  used_restart_global = nullptr;
  nfix_restart_peratom = 0;
  id_restart_peratom = style_restart_peratom = nullptr;
  index_restart_peratom = used_restart_peratom = nullptr;
  create_factories();
}
void _noopt Modify::create_factories() {
  fix_map = new FixCreatorMap();
#define FIX_CLASS
#define FixStyle(key, Class) (*fix_map)[#key] = &style_creator<Fix, Class>;
#include "style_fix.h"
#undef FixStyle
#undef FIX_CLASS
}
Modify::~Modify() {
  while (nfix)
    delete_fix(0);
  memory->sfree(fix);
  memory->destroy(fmask);
  delete[] list_initial_integrate;
  delete[] list_post_integrate;
  delete[] list_pre_exchange;
  delete[] list_pre_neighbor;
  delete[] list_post_neighbor;
  delete[] list_pre_force;
  delete[] list_pre_reverse;
  delete[] list_post_force;
  delete[] list_post_force_group;
  delete[] list_final_integrate;
  delete[] list_end_of_step;
  delete[] list_energy_couple;
  delete[] list_energy_global;
  delete[] list_energy_atom;
  delete[] list_initial_integrate_respa;
  delete[] list_post_integrate_respa;
  delete[] list_pre_force_respa;
  delete[] list_post_force_respa;
  delete[] list_final_integrate_respa;
  delete[] list_min_pre_exchange;
  delete[] list_min_pre_neighbor;
  delete[] list_min_post_neighbor;
  delete[] list_min_pre_force;
  delete[] list_min_pre_reverse;
  delete[] list_min_post_force;
  delete[] list_min_energy;
  delete[] end_of_step_every;
  delete[] list_timeflag;
  delete fix_map;
}
void Modify::init() {
  int i, j;
  for (i = 0; i < nfix; i++)
    fix[i]->init();
  restart_pbc_any = 0;
  for (i = 0; i < nfix; i++)
    if (fix[i]->restart_pbc)
      restart_pbc_any = 1;
  list_init(INITIAL_INTEGRATE, n_initial_integrate, list_initial_integrate);
  list_init(POST_INTEGRATE, n_post_integrate, list_post_integrate);
  list_init(PRE_EXCHANGE, n_pre_exchange, list_pre_exchange);
  list_init(PRE_NEIGHBOR, n_pre_neighbor, list_pre_neighbor);
  list_init(POST_NEIGHBOR, n_post_neighbor, list_post_neighbor);
  list_init(PRE_FORCE, n_pre_force, list_pre_force);
  list_init(PRE_REVERSE, n_pre_reverse, list_pre_reverse);
  list_init(POST_FORCE, n_post_force, list_post_force);
  list_init_post_force_group(n_post_force_group, list_post_force_group);
  list_init(FINAL_INTEGRATE, n_final_integrate, list_final_integrate);
  list_init_end_of_step(END_OF_STEP, n_end_of_step, list_end_of_step);
  list_init_energy_couple(n_energy_couple, list_energy_couple);
  list_init_energy_global(n_energy_global, list_energy_global);
  list_init_energy_atom(n_energy_atom, list_energy_atom);
  list_init(INITIAL_INTEGRATE_RESPA, n_initial_integrate_respa,
            list_initial_integrate_respa);
  list_init(POST_INTEGRATE_RESPA, n_post_integrate_respa,
            list_post_integrate_respa);
  list_init(POST_FORCE_RESPA, n_post_force_respa, list_post_force_respa);
  list_init(PRE_FORCE_RESPA, n_pre_force_respa, list_pre_force_respa);
  list_init(FINAL_INTEGRATE_RESPA, n_final_integrate_respa,
            list_final_integrate_respa);
  list_init(MIN_PRE_EXCHANGE, n_min_pre_exchange, list_min_pre_exchange);
  list_init(MIN_PRE_NEIGHBOR, n_min_pre_neighbor, list_min_pre_neighbor);
  list_init(MIN_POST_NEIGHBOR, n_min_post_neighbor, list_min_post_neighbor);
  list_init(MIN_PRE_FORCE, n_min_pre_force, list_min_pre_force);
  list_init(MIN_PRE_REVERSE, n_min_pre_reverse, list_min_pre_reverse);
  list_init(MIN_POST_FORCE, n_min_post_force, list_min_post_force);
  list_init(MIN_ENERGY, n_min_energy, list_min_energy);
  n_post_force_any = n_post_force + n_post_force_group;
  n_post_force_respa_any = n_post_force_respa + n_post_force_group;
  for (i = 0; i < nfix; i++)
    if (!fix[i]->dynamic_group_allow && group->dynamic[fix[i]->igroup])
      error->all(FLERR, "Fix {} does not allow use with a dynamic group",
                 fix[i]->style);
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *flag = new int[nlocal];
  for (i = 0; i < nlocal; i++)
    flag[i] = 0;
  int groupbit;
  for (i = 0; i < nfix; i++) {
    if (fix[i]->time_integrate == 0)
      continue;
    groupbit = fix[i]->groupbit;
    for (j = 0; j < nlocal; j++)
      if (mask[j] & groupbit)
        flag[j]++;
  }
  int check = 0;
  for (i = 0; i < nlocal; i++)
    if (flag[i] > 1)
      check = 1;
  delete[] flag;
  int checkall;
  MPI_Allreduce(&check, &checkall, 1, MPI_INT, MPI_SUM, world);
  if (comm->me == 0 && checkall)
    error->warning(FLERR,
                   "One or more atoms are time integrated more than once");
}
void Modify::setup(int vflag) {
  for (int i = 0; i < nfix; i++)
    if (strcmp(fix[i]->style, "GROUP") == 0)
      fix[i]->setup(vflag);
  if (update->whichflag == 1)
    for (int i = 0; i < nfix; i++)
      fix[i]->setup(vflag);
  else if (update->whichflag == 2)
    for (int i = 0; i < nfix; i++)
      fix[i]->min_setup(vflag);
}
void Modify::setup_pre_exchange() {
  if (update->whichflag <= 1)
    for (int i = 0; i < n_pre_exchange; i++)
      fix[list_pre_exchange[i]]->setup_pre_exchange();
  else if (update->whichflag == 2)
    for (int i = 0; i < n_min_pre_exchange; i++)
      fix[list_min_pre_exchange[i]]->setup_pre_exchange();
}
void Modify::setup_pre_neighbor() {
  if (update->whichflag == 1)
    for (int i = 0; i < n_pre_neighbor; i++)
      fix[list_pre_neighbor[i]]->setup_pre_neighbor();
  else if (update->whichflag == 2)
    for (int i = 0; i < n_min_pre_neighbor; i++)
      fix[list_min_pre_neighbor[i]]->setup_pre_neighbor();
}
void Modify::setup_post_neighbor() {
  if (update->whichflag == 1)
    for (int i = 0; i < n_post_neighbor; i++)
      fix[list_post_neighbor[i]]->setup_post_neighbor();
  else if (update->whichflag == 2)
    for (int i = 0; i < n_min_post_neighbor; i++)
      fix[list_min_post_neighbor[i]]->setup_post_neighbor();
}
void Modify::setup_pre_force(int vflag) {
  if (update->whichflag == 1)
    for (int i = 0; i < n_pre_force; i++)
      fix[list_pre_force[i]]->setup_pre_force(vflag);
  else if (update->whichflag == 2)
    for (int i = 0; i < n_min_pre_force; i++)
      fix[list_min_pre_force[i]]->setup_pre_force(vflag);
}
void Modify::setup_pre_reverse(int eflag, int vflag) {
  if (update->whichflag == 1)
    for (int i = 0; i < n_pre_reverse; i++)
      fix[list_pre_reverse[i]]->setup_pre_reverse(eflag, vflag);
  else if (update->whichflag == 2)
    for (int i = 0; i < n_min_pre_reverse; i++)
      fix[list_min_pre_reverse[i]]->setup_pre_reverse(eflag, vflag);
}
void Modify::initial_integrate(int vflag) {
  for (int i = 0; i < n_initial_integrate; i++)
    fix[list_initial_integrate[i]]->initial_integrate(vflag);
}
void Modify::post_integrate() {
  for (int i = 0; i < n_post_integrate; i++)
    fix[list_post_integrate[i]]->post_integrate();
}
void Modify::pre_exchange() {
  for (int i = 0; i < n_pre_exchange; i++)
    fix[list_pre_exchange[i]]->pre_exchange();
}
void Modify::pre_neighbor() {
  for (int i = 0; i < n_pre_neighbor; i++)
    fix[list_pre_neighbor[i]]->pre_neighbor();
}
void Modify::post_neighbor() {
  for (int i = 0; i < n_post_neighbor; i++)
    fix[list_post_neighbor[i]]->post_neighbor();
}
void Modify::pre_force(int vflag) {
  for (int i = 0; i < n_pre_force; i++)
    fix[list_pre_force[i]]->pre_force(vflag);
}
void Modify::pre_reverse(int eflag, int vflag) {
  for (int i = 0; i < n_pre_reverse; i++)
    fix[list_pre_reverse[i]]->pre_reverse(eflag, vflag);
}
void Modify::post_force(int vflag) {
  if (n_post_force_group) {
    for (int i = 0; i < n_post_force_group; i++)
      fix[list_post_force_group[i]]->post_force(vflag);
  }
  if (n_post_force) {
    for (int i = 0; i < n_post_force; i++)
      fix[list_post_force[i]]->post_force(vflag);
  }
}
void Modify::final_integrate() {
  for (int i = 0; i < n_final_integrate; i++)
    fix[list_final_integrate[i]]->final_integrate();
}
void Modify::end_of_step() {
  for (int i = 0; i < n_end_of_step; i++)
    if (update->ntimestep % end_of_step_every[i] == 0)
      fix[list_end_of_step[i]]->end_of_step();
}
double Modify::energy_global() {
  double energy = 0.0;
  for (int i = 0; i < n_energy_global; i++)
    energy += fix[list_energy_global[i]]->compute_scalar();
  return energy;
}
void Modify::energy_atom(int nlocal, double *energy) {
  int i, j;
  double *eatom;
  for (i = 0; i < n_energy_atom; i++) {
    eatom = fix[list_energy_atom[i]]->eatom;
    if (!eatom)
      continue;
    for (j = 0; j < nlocal; j++)
      energy[j] += eatom[j];
  }
}
void Modify::post_run() {
  for (int i = 0; i < nfix; i++)
    fix[i]->post_run();
  n_timeflag = -1;
}
void Modify::setup_pre_force_respa(int vflag, int ilevel) {
  for (int i = 0; i < n_pre_force_respa; i++)
    fix[list_pre_force_respa[i]]->setup_pre_force_respa(vflag, ilevel);
}
void Modify::initial_integrate_respa(int vflag, int ilevel, int iloop) {
  for (int i = 0; i < n_initial_integrate_respa; i++)
    fix[list_initial_integrate_respa[i]]->initial_integrate_respa(vflag, ilevel,
                                                                  iloop);
}
void Modify::post_integrate_respa(int ilevel, int iloop) {
  for (int i = 0; i < n_post_integrate_respa; i++)
    fix[list_post_integrate_respa[i]]->post_integrate_respa(ilevel, iloop);
}
void Modify::pre_force_respa(int vflag, int ilevel, int iloop) {
  for (int i = 0; i < n_pre_force_respa; i++)
    fix[list_pre_force_respa[i]]->pre_force_respa(vflag, ilevel, iloop);
}
void Modify::post_force_respa(int vflag, int ilevel, int iloop) {
  if (n_post_force_group) {
    for (int i = 0; i < n_post_force_group; i++)
      fix[list_post_force_group[i]]->post_force_respa(vflag, ilevel, iloop);
  }
  if (n_post_force_respa) {
    for (int i = 0; i < n_post_force_respa; i++)
      fix[list_post_force_respa[i]]->post_force_respa(vflag, ilevel, iloop);
  }
}
void Modify::final_integrate_respa(int ilevel, int iloop) {
  for (int i = 0; i < n_final_integrate_respa; i++)
    fix[list_final_integrate_respa[i]]->final_integrate_respa(ilevel, iloop);
}
void Modify::min_pre_exchange() {
  for (int i = 0; i < n_min_pre_exchange; i++)
    fix[list_min_pre_exchange[i]]->min_pre_exchange();
}
void Modify::min_pre_neighbor() {
  for (int i = 0; i < n_min_pre_neighbor; i++)
    fix[list_min_pre_neighbor[i]]->min_pre_neighbor();
}
void Modify::min_post_neighbor() {
  for (int i = 0; i < n_min_post_neighbor; i++)
    fix[list_min_post_neighbor[i]]->min_post_neighbor();
}
void Modify::min_pre_force(int vflag) {
  for (int i = 0; i < n_min_pre_force; i++)
    fix[list_min_pre_force[i]]->min_pre_force(vflag);
}
void Modify::min_pre_reverse(int eflag, int vflag) {
  for (int i = 0; i < n_min_pre_reverse; i++)
    fix[list_min_pre_reverse[i]]->min_pre_reverse(eflag, vflag);
}
void Modify::min_post_force(int vflag) {
  for (int i = 0; i < n_min_post_force; i++)
    fix[list_min_post_force[i]]->min_post_force(vflag);
}
double Modify::min_energy(double *fextra) {
  int ifix, index;
  index = 0;
  double eng = 0.0;
  for (int i = 0; i < n_min_energy; i++) {
    ifix = list_min_energy[i];
    eng += fix[ifix]->min_energy(&fextra[index]);
    index += fix[ifix]->min_dof();
  }
  return eng;
}
void Modify::min_store() {
  for (int i = 0; i < n_min_energy; i++)
    fix[list_min_energy[i]]->min_store();
}
void Modify::min_clearstore() {
  for (int i = 0; i < n_min_energy; i++)
    fix[list_min_energy[i]]->min_clearstore();
}
void Modify::min_pushstore() {
  for (int i = 0; i < n_min_energy; i++)
    fix[list_min_energy[i]]->min_pushstore();
}
void Modify::min_popstore() {
  for (int i = 0; i < n_min_energy; i++)
    fix[list_min_energy[i]]->min_popstore();
}
void Modify::min_step(double alpha, double *hextra) {
  int ifix, index;
  index = 0;
  for (int i = 0; i < n_min_energy; i++) {
    ifix = list_min_energy[i];
    fix[ifix]->min_step(alpha, &hextra[index]);
    index += fix[ifix]->min_dof();
  }
}
double Modify::max_alpha(double *hextra) {
  int ifix, index;
  double alpha = BIG;
  index = 0;
  for (int i = 0; i < n_min_energy; i++) {
    ifix = list_min_energy[i];
    double alpha_one = fix[ifix]->max_alpha(&hextra[index]);
    alpha = MIN(alpha, alpha_one);
    index += fix[ifix]->min_dof();
  }
  return alpha;
}
int Modify::min_dof() {
  int ndof = 0;
  for (int i = 0; i < n_min_energy; i++)
    ndof += fix[list_min_energy[i]]->min_dof();
  return ndof;
}
int Modify::min_reset_ref() {
  int itmp, itmpall;
  itmpall = 0;
  for (int i = 0; i < n_min_energy; i++) {
    itmp = fix[list_min_energy[i]]->min_reset_ref();
    if (itmp)
      itmpall = 1;
  }
  return itmpall;
}
void Modify::reset_grid() {
  for (int i = 0; i < nfix; i++)
    if (fix[i]->pergrid_flag)
      fix[i]->reset_grid();
}
Fix *Modify::add_fix(int narg, char **arg, int trysuffix) {
  if (narg < 3)
    utils::missing_cmd_args(FLERR, "fix", error);
  const char *exceptions[] = {"GPU",
                              "OMP",
                              "INTEL",
                              "property/atom",
                              "cmap",
                              "cmap3",
                              "rx",
                              "deprecated",
                              "STORE/KIM",
                              "amoeba/pitorsion",
                              "amoeba/bitorsion",
                              nullptr};
  if (domain->box_exist == 0) {
    int m;
    for (m = 0; exceptions[m] != nullptr; m++)
      if (strcmp(arg[2], exceptions[m]) == 0)
        break;
    if (exceptions[m] == nullptr)
      error->all(FLERR, "Fix command before simulation box is defined");
  }
  int igroup = group->find(arg[1]);
  if (igroup == -1)
    error->all(FLERR, "Could not find fix group ID {}", arg[1]);
  int ifix, newflag;
  for (ifix = 0; ifix < nfix; ifix++)
    if (strcmp(arg[0], fix[ifix]->id) == 0)
      break;
  if (ifix < nfix) {
    newflag = 0;
    int match = 0;
    if (strcmp(arg[2], fix[ifix]->style) == 0)
      match = 1;
    if (!match)
      error->all(FLERR, "Replacing a fix, but new style != old style");
    if (fix[ifix]->igroup != igroup && comm->me == 0)
      error->warning(FLERR, "Replacing a fix, but new group != old group");
    delete fix[ifix];
    fix[ifix] = nullptr;
  } else {
    newflag = 1;
    if (nfix == maxfix) {
      maxfix += DELTA;
      fix = (Fix **)memory->srealloc(fix, maxfix * sizeof(Fix *), "modify:fix");
      memory->grow(fmask, maxfix, "modify:fmask");
    }
  }
  fix[ifix] = nullptr;
  if ((fix[ifix] == nullptr) && (fix_map->find(arg[2]) != fix_map->end())) {
    FixCreator &fix_creator = (*fix_map)[arg[2]];
    fix[ifix] = fix_creator(lmp, narg, arg);
  }
  if (newflag) {
    nfix++;
    fix_list = std::vector<Fix *>(fix, fix + nfix);
  }
  fix[ifix]->post_constructor();
  for (int i = 0; i < nfix_restart_global; i++)
    if ((strcmp(id_restart_global[i], fix[ifix]->id) == 0) &&
        (utils::strip_style_suffix(fix[ifix]->style, lmp) ==
         style_restart_global[i])) {
      fix[ifix]->restart(state_restart_global[i]);
      used_restart_global[i] = 1;
      fix[ifix]->restart_reset = 1;
      if (comm->me == 0)
        utils::logmesg(lmp,
                       "Resetting global fix info from restart file:\n"
                       "  fix style: {}, fix ID: {}\n",
                       fix[ifix]->style, fix[ifix]->id);
    }
  for (int i = 0; i < nfix_restart_peratom; i++)
    if (strcmp(id_restart_peratom[i], fix[ifix]->id) == 0 &&
        strcmp(style_restart_peratom[i], fix[ifix]->style) == 0) {
      used_restart_peratom[i] = 1;
      for (int j = 0; j < atom->nlocal; j++)
        fix[ifix]->unpack_restart(j, index_restart_peratom[i]);
      fix[ifix]->restart_reset = 1;
      if (comm->me == 0)
        utils::logmesg(lmp,
                       "Resetting peratom fix info from restart file:\n"
                       "  fix style: {}, fix ID: {}\n",
                       fix[ifix]->style, fix[ifix]->id);
    }
  fmask[ifix] = fix[ifix]->setmask();
  return fix[ifix];
}
Fix *Modify::add_fix(const std::string &fixcmd, int trysuffix) {
  auto args = utils::split_words(fixcmd);
  std::vector<char *> newarg(args.size());
  int i = 0;
  for (const auto &arg : args) {
    newarg[i++] = (char *)arg.c_str();
  }
  return add_fix(args.size(), newarg.data(), trysuffix);
}
void Modify::modify_fix(int narg, char **arg) {
  if (narg < 2)
    utils::missing_cmd_args(FLERR, "fix_modify", error);
  auto ifix = get_fix_by_id(arg[0]);
  if (!ifix)
    error->all(FLERR, "Could not find fix_modify ID {}", arg[0]);
  ifix->modify_params(narg - 1, &arg[1]);
}
void Modify::delete_fix(const std::string &id) {
  int ifix = find_fix(id);
  if (ifix < 0)
    error->all(FLERR, "Could not find fix ID {} to delete", id);
  delete_fix(ifix);
}
void Modify::delete_fix(int ifix) {
  if ((ifix < 0) || (ifix >= nfix))
    return;
  delete fix[ifix];
  for (int i = ifix + 1; i < nfix; i++)
    fix[i - 1] = fix[i];
  for (int i = ifix + 1; i < nfix; i++)
    fmask[i - 1] = fmask[i];
  nfix--;
  fix_list = std::vector<Fix *>(fix, fix + nfix);
}
int Modify::find_fix(const std::string &id) {
  if (id.empty())
    return -1;
  for (int ifix = 0; ifix < nfix; ifix++)
    if (fix[ifix] && (id == fix[ifix]->id))
      return ifix;
  return -1;
}
Fix *Modify::get_fix_by_id(const std::string &id) const {
  if (id.empty())
    return nullptr;
  for (int ifix = 0; ifix < nfix; ifix++)
    if (fix[ifix] && (id == fix[ifix]->id))
      return fix[ifix];
  return nullptr;
}
const std::vector<Fix *> &Modify::get_fix_list() {
  fix_list = std::vector<Fix *>(fix, fix + nfix);
  return fix_list;
}
void Modify::list_init(int mask, int &n, int *&list) {
  delete[] list;
  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fmask[i] & mask)
      n++;
  list = new int[n];
  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fmask[i] & mask)
      list[n++] = i;
}
void Modify::list_init_end_of_step(int mask, int &n, int *&list) {
  delete[] list;
  delete[] end_of_step_every;
  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fmask[i] & mask)
      n++;
  list = new int[n];
  end_of_step_every = new int[n];
  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fmask[i] & mask) {
      list[n] = i;
      end_of_step_every[n++] = fix[i]->nevery;
    }
}
void Modify::list_init_energy_couple(int &n, int *&list) {
  delete[] list;
  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fix[i]->ecouple_flag)
      n++;
  list = new int[n];
  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fix[i]->ecouple_flag)
      list[n++] = i;
}
void Modify::list_init_energy_global(int &n, int *&list) {
  delete[] list;
  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fix[i]->energy_global_flag && fix[i]->thermo_energy)
      n++;
  list = new int[n];
  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fix[i]->energy_global_flag && fix[i]->thermo_energy)
      list[n++] = i;
}
void Modify::list_init_energy_atom(int &n, int *&list) {
  delete[] list;
  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fix[i]->energy_peratom_flag && fix[i]->thermo_energy)
      n++;
  list = new int[n];
  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fix[i]->energy_peratom_flag && fix[i]->thermo_energy)
      list[n++] = i;
}
void Modify::list_init_post_force_group(int &n, int *&list) {
  delete[] list;
  n = 0;
  for (int i = 0; i < nfix; i++)
    if (strcmp(fix[i]->style, "GROUP") == 0)
      n++;
  list = new int[n];
  n = 0;
  for (int i = 0; i < nfix; i++)
    if (strcmp(fix[i]->style, "GROUP") == 0)
      list[n++] = i;
}
