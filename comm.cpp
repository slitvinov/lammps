#include <map>
#include <set>
#include <vector>
#include <unordered_set>
#include <string>
#include <cmath>
#include <mpi.h>
#include "lmptype.h"
#include "utils.h"
#include "lammps.h"
#include "pointers.h"
#include "comm.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "pair.h"
#include "procmap.h"
#include "universe.h"
#include "update.h"
#include <cstring>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace LAMMPS_NS;
#define BUFEXTRA 1024
enum { ONELEVEL, TWOLEVEL, NUMA, CUSTOM };
enum { CART, CARTREORDER, XYZ };
Comm::Comm(LAMMPS *lmp) : Pointers(lmp) {
  MPI_Comm_rank(world, &me);
  MPI_Comm_size(world, &nprocs);
  mode = 0;
  bordergroup = 0;
  cutghostuser = 0.0;
  cutusermulti = nullptr;
  cutusermultiold = nullptr;
  ncollections = 0;
  ncollections_cutoff = 0;
  ghost_velocity = 0;
  user_procgrid[0] = user_procgrid[1] = user_procgrid[2] = 0;
  coregrid[0] = coregrid[1] = coregrid[2] = 1;
  gridflag = ONELEVEL;
  mapflag = CART;
  customfile = nullptr;
  outfile = nullptr;
  recv_from_partition = send_to_partition = -1;
  otherflag = 0;
  maxexchange = maxexchange_atom = maxexchange_fix = 0;
  maxexchange_fix_dynamic = 0;
  bufextra = BUFEXTRA;
  grid2proc = nullptr;
  xsplit = ysplit = zsplit = nullptr;
  rcbnew = 0;
  multi_reduce = 0;
  nthreads = 1;
#ifdef _OPENMP
  if (getenv("OMP_NUM_THREADS") == nullptr) {
    nthreads = 1;
    if (me == 0)
      error->message(FLERR, "OMP_NUM_THREADS environment is not set. "
                            "Defaulting to 1 thread.");
  } else {
    nthreads = omp_get_max_threads();
  }
  MPI_Bcast(&nthreads, 1, MPI_INT, 0, world);
  omp_set_num_threads(nthreads);
  if (me == 0)
    utils::logmesg(lmp, "  using {} OpenMP thread(s) per MPI task\n", nthreads);
#endif
}
void Comm::init() {
  triclinic = domain->triclinic;
  map_style = atom->map_style;
  domain->subbox_too_small_check(neighbor->skin);
  comm_x_only = atom->avec->comm_x_only;
  comm_f_only = atom->avec->comm_f_only;
  if (ghost_velocity)
    comm_x_only = 0;
  size_forward = atom->avec->size_forward;
  size_reverse = atom->avec->size_reverse;
  size_border = atom->avec->size_border;
  if (ghost_velocity)
    size_forward += atom->avec->size_velocity;
  if (ghost_velocity)
    size_border += atom->avec->size_velocity;
  const auto &fix_list = modify->get_fix_list();
  for (const auto &fix : fix_list)
    size_border += fix->comm_border;
  maxforward = MAX(size_forward, size_border);
  maxreverse = size_reverse;
  if (force->pair)
    maxforward = MAX(maxforward, force->pair->comm_forward);
  if (force->pair)
    maxreverse = MAX(maxreverse, force->pair->comm_reverse);
  for (const auto &fix : fix_list) {
    maxforward = MAX(maxforward, fix->comm_forward);
    maxreverse = MAX(maxreverse, fix->comm_reverse);
  }
  if (force->newton == 0)
    maxreverse = 0;
  if (force->pair)
    maxreverse = MAX(maxreverse, force->pair->comm_reverse_off);
  maxexchange_atom = atom->avec->maxexchange;
  maxexchange_fix_dynamic = 0;
  for (const auto &fix : fix_list)
    if (fix->maxexchange_dynamic)
      maxexchange_fix_dynamic = 1;
}
void Comm::init_exchange() {
  maxexchange_fix = 0;
  for (const auto &fix : modify->get_fix_list())
    maxexchange_fix += fix->maxexchange;
  maxexchange = maxexchange_atom + maxexchange_fix;
  bufextra = maxexchange + BUFEXTRA;
}
void Comm::modify_params(int narg, char **arg) {
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "vel") == 0) {
      ghost_velocity = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    }
  }
}
void Comm::set_proc_grid(int outflag) {
  if (recv_from_partition >= 0) {
    if (me == 0) {
      MPI_Recv(other_procgrid, 3, MPI_INT,
               universe->root_proc[recv_from_partition], 0, universe->uworld,
               MPI_STATUS_IGNORE);
      MPI_Recv(other_coregrid, 3, MPI_INT,
               universe->root_proc[recv_from_partition], 0, universe->uworld,
               MPI_STATUS_IGNORE);
    }
    MPI_Bcast(other_procgrid, 3, MPI_INT, 0, world);
    MPI_Bcast(other_coregrid, 3, MPI_INT, 0, world);
  }
  auto pmap = new ProcMap(lmp);
  if (gridflag == ONELEVEL) {
    pmap->onelevel_grid(nprocs, user_procgrid, procgrid, otherflag, other_style,
                        other_procgrid, other_coregrid);
  }
  if (grid2proc)
    memory->destroy(grid2proc);
  memory->create(grid2proc, procgrid[0], procgrid[1], procgrid[2],
                 "comm:grid2proc");
  if (gridflag == ONELEVEL) {
    if (mapflag == CART)
      pmap->cart_map(0, procgrid, myloc, procneigh, grid2proc);
    else if (mapflag == CARTREORDER)
      pmap->cart_map(1, procgrid, myloc, procneigh, grid2proc);
  }
  delete pmap;
  memory->destroy(xsplit);
  memory->destroy(ysplit);
  memory->destroy(zsplit);
  memory->create(xsplit, procgrid[0] + 1, "comm:xsplit");
  memory->create(ysplit, procgrid[1] + 1, "comm:ysplit");
  memory->create(zsplit, procgrid[2] + 1, "comm:zsplit");
  for (int i = 0; i < procgrid[0]; i++)
    xsplit[i] = i * 1.0 / procgrid[0];
  for (int i = 0; i < procgrid[1]; i++)
    ysplit[i] = i * 1.0 / procgrid[1];
  for (int i = 0; i < procgrid[2]; i++)
    zsplit[i] = i * 1.0 / procgrid[2];
  xsplit[procgrid[0]] = ysplit[procgrid[1]] = zsplit[procgrid[2]] = 1.0;
  if (send_to_partition >= 0) {
    if (me == 0) {
      MPI_Send(procgrid, 3, MPI_INT, universe->root_proc[send_to_partition], 0,
               universe->uworld);
      MPI_Send(coregrid, 3, MPI_INT, universe->root_proc[send_to_partition], 0,
               universe->uworld);
    }
  }
}
double Comm::get_comm_cutoff() {
  double maxcommcutoff, maxbondcutoff = 0.0;
  maxcommcutoff = MAX(cutghostuser, neighbor->cutneighmax);
  if (!force->pair && (cutghostuser == 0.0)) {
    maxcommcutoff = MAX(maxcommcutoff, maxbondcutoff);
  }
  if (neighbor->interval_collection_flag) {
    for (int i = 0; i < neighbor->ncollections; i++) {
      maxcommcutoff = MAX(maxcommcutoff, neighbor->collection2cut[i]);
    }
  }
  return maxcommcutoff;
}
