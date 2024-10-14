#include "comm.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "error.h"
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
Comm::~Comm() {
  memory->destroy(grid2proc);
  memory->destroy(xsplit);
  memory->destroy(ysplit);
  memory->destroy(zsplit);
  memory->destroy(cutusermulti);
  memory->destroy(cutusermultiold);
  delete[] customfile;
  delete[] outfile;
}
void Comm::copy_arrays(Comm *oldcomm) {
  if (oldcomm->grid2proc) {
    memory->create(grid2proc, procgrid[0], procgrid[1], procgrid[2],
                   "comm:grid2proc");
    memcpy(&grid2proc[0][0][0], &oldcomm->grid2proc[0][0][0],
           (procgrid[0] * procgrid[1] * procgrid[2]) * sizeof(int));
    memory->create(xsplit, procgrid[0] + 1, "comm:xsplit");
    memory->create(ysplit, procgrid[1] + 1, "comm:ysplit");
    memory->create(zsplit, procgrid[2] + 1, "comm:zsplit");
    memcpy(xsplit, oldcomm->xsplit, (procgrid[0] + 1) * sizeof(double));
    memcpy(ysplit, oldcomm->ysplit, (procgrid[1] + 1) * sizeof(double));
    memcpy(zsplit, oldcomm->zsplit, (procgrid[2] + 1) * sizeof(double));
  }
  ncollections = oldcomm->ncollections;
  ncollections_cutoff = oldcomm->ncollections_cutoff;
  if (oldcomm->cutusermulti) {
    memory->create(cutusermulti, ncollections_cutoff, "comm:cutusermulti");
    memcpy(cutusermulti, oldcomm->cutusermulti, ncollections_cutoff);
  }
  if (oldcomm->cutusermultiold) {
    memory->create(cutusermultiold, atom->ntypes + 1, "comm:cutusermultiold");
    memcpy(cutusermultiold, oldcomm->cutusermultiold, atom->ntypes + 1);
  }
  if (customfile)
    customfile = utils::strdup(oldcomm->customfile);
  if (outfile)
    outfile = utils::strdup(oldcomm->outfile);
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
  if ((mode == Comm::MULTI) && (neighbor->style != Neighbor::MULTI))
    error->all(FLERR,
               "Cannot use comm mode multi without multi-style neighbor lists");
}
void Comm::init_exchange() {
  maxexchange_fix = 0;
  for (const auto &fix : modify->get_fix_list())
    maxexchange_fix += fix->maxexchange;
  maxexchange = maxexchange_atom + maxexchange_fix;
  bufextra = maxexchange + BUFEXTRA;
}
void Comm::modify_params(int narg, char **arg) {
  if (narg < 1)
    utils::missing_cmd_args(FLERR, "comm_modify", error);
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "vel") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, "comm_modify vel", error);
      ghost_velocity = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else
      error->all(FLERR, "Unknown comm_modify keyword: {}", arg[iarg]);
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
  } else if (gridflag == TWOLEVEL) {
    pmap->twolevel_grid(nprocs, user_procgrid, procgrid, ncores, user_coregrid,
                        coregrid, otherflag, other_style, other_procgrid,
                        other_coregrid);
  } else if (gridflag == NUMA) {
    pmap->numa_grid(nprocs, user_procgrid, procgrid, coregrid);
  } else if (gridflag == CUSTOM) {
    pmap->custom_grid(customfile, nprocs, user_procgrid, procgrid);
  }
  if (procgrid[0] * procgrid[1] * procgrid[2] != nprocs)
    error->all(FLERR, "Bad grid of processors");
  if (domain->dimension == 2 && procgrid[2] != 1)
    error->all(FLERR, "Processor count in z must be 1 for 2d simulation");
  if (grid2proc)
    memory->destroy(grid2proc);
  memory->create(grid2proc, procgrid[0], procgrid[1], procgrid[2],
                 "comm:grid2proc");
  if (gridflag == ONELEVEL) {
    if (mapflag == CART)
      pmap->cart_map(0, procgrid, myloc, procneigh, grid2proc);
    else if (mapflag == CARTREORDER)
      pmap->cart_map(1, procgrid, myloc, procneigh, grid2proc);
    else if (mapflag == XYZ)
      pmap->xyz_map(xyz, procgrid, myloc, procneigh, grid2proc);
  } else if (gridflag == TWOLEVEL) {
    if (mapflag == CART)
      pmap->cart_map(0, procgrid, ncores, coregrid, myloc, procneigh,
                     grid2proc);
    else if (mapflag == CARTREORDER)
      pmap->cart_map(1, procgrid, ncores, coregrid, myloc, procneigh,
                     grid2proc);
    else if (mapflag == XYZ)
      pmap->xyz_map(xyz, procgrid, ncores, coregrid, myloc, procneigh,
                    grid2proc);
  } else if (gridflag == NUMA) {
    pmap->numa_map(0, coregrid, myloc, procneigh, grid2proc);
  } else if (gridflag == CUSTOM) {
    pmap->custom_map(procgrid, myloc, procneigh, grid2proc);
  }
  if (outflag && me == 0) {
    auto mesg = fmt::format("  {} by {} by {} MPI processor grid\n",
                            procgrid[0], procgrid[1], procgrid[2]);
    if (gridflag == NUMA || gridflag == TWOLEVEL)
      mesg += fmt::format("  {} by {} by {} core grid within node\n",
                          coregrid[0], coregrid[1], coregrid[2]);
    utils::logmesg(lmp, mesg);
  }
  if (outfile)
    pmap->output(outfile, procgrid, grid2proc);
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
  if (domain->triclinic)
    domain->set_lamda_box();
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
  } else {
    if ((me == 0) && (maxbondcutoff > maxcommcutoff))
      error->warning(FLERR,
                     "Communication cutoff {} is shorter than a bond "
                     "length based estimate of {}. This may lead to errors.",
                     maxcommcutoff, maxbondcutoff);
  }
  if ((me == 0) && (update->setupflag == 1)) {
    if ((cutghostuser > 0.0) && (maxcommcutoff > cutghostuser))
      error->warning(FLERR, "Communication cutoff adjusted to {}",
                     maxcommcutoff);
  }
  if (neighbor->interval_collection_flag) {
    for (int i = 0; i < neighbor->ncollections; i++) {
      maxcommcutoff = MAX(maxcommcutoff, neighbor->collection2cut[i]);
    }
  }
  return maxcommcutoff;
}
int Comm::coord2proc(double *x, int &igx, int &igy, int &igz) {
  double *prd = domain->prd;
  double *boxlo = domain->boxlo;
  triclinic = domain->triclinic;
  if (layout == Comm::LAYOUT_UNIFORM) {
    if (triclinic == 0) {
      igx = static_cast<int>(procgrid[0] * (x[0] - boxlo[0]) / prd[0]);
      igy = static_cast<int>(procgrid[1] * (x[1] - boxlo[1]) / prd[1]);
      igz = static_cast<int>(procgrid[2] * (x[2] - boxlo[2]) / prd[2]);
    } else {
      igx = static_cast<int>(procgrid[0] * x[0]);
      igy = static_cast<int>(procgrid[1] * x[1]);
      igz = static_cast<int>(procgrid[2] * x[2]);
    }
  } else if (layout == Comm::LAYOUT_NONUNIFORM) {
    if (triclinic == 0) {
      igx =
          utils::binary_search((x[0] - boxlo[0]) / prd[0], procgrid[0], xsplit);
      igy =
          utils::binary_search((x[1] - boxlo[1]) / prd[1], procgrid[1], ysplit);
      igz =
          utils::binary_search((x[2] - boxlo[2]) / prd[2], procgrid[2], zsplit);
    } else {
      igx = utils::binary_search(x[0], procgrid[0], xsplit);
      igy = utils::binary_search(x[1], procgrid[1], ysplit);
      igz = utils::binary_search(x[2], procgrid[2], zsplit);
    }
  }
  if (igx < 0)
    igx = 0;
  if (igx >= procgrid[0])
    igx = procgrid[0] - 1;
  if (igy < 0)
    igy = 0;
  if (igy >= procgrid[1])
    igy = procgrid[1] - 1;
  if (igz < 0)
    igz = 0;
  if (igz >= procgrid[2])
    igz = procgrid[2] - 1;
  return grid2proc[igx][igy][igz];
}
void Comm::ring(int n, int nper, void *inbuf, int messtag,
                void (*callback)(int, char *, void *), void *outbuf, void *ptr,
                int self) {
  MPI_Request request;
  MPI_Status status;
  int nbytes = n * nper;
  int maxbytes;
  MPI_Allreduce(&nbytes, &maxbytes, 1, MPI_INT, MPI_MAX, world);
  if (maxbytes == 0)
    return;
  if ((nbytes > 0) && inbuf == nullptr)
    error->one(FLERR, "Cannot put data on ring from NULL pointer");
  char *buf, *bufcopy;
  memory->create(buf, maxbytes, "comm:buf");
  memory->create(bufcopy, maxbytes, "comm:bufcopy");
  if (nbytes && inbuf)
    memcpy(buf, inbuf, nbytes);
  int next = me + 1;
  int prev = me - 1;
  if (next == nprocs)
    next = 0;
  if (prev < 0)
    prev = nprocs - 1;
  for (int loop = 0; loop < nprocs; loop++) {
    if (me != next) {
      MPI_Irecv(bufcopy, maxbytes, MPI_CHAR, prev, messtag, world, &request);
      MPI_Send(buf, nbytes, MPI_CHAR, next, messtag, world);
      MPI_Wait(&request, &status);
      MPI_Get_count(&status, MPI_CHAR, &nbytes);
      if (nbytes)
        memcpy(buf, bufcopy, nbytes);
    }
    if (self || loop < nprocs - 1)
      callback(nbytes / nper, buf, ptr);
  }
  if (nbytes && outbuf)
    memcpy(outbuf, buf, nbytes);
  memory->destroy(buf);
  memory->destroy(bufcopy);
}
