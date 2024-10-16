#include "comm_brick.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "memory.h"
#include "neighbor.h"
#include "pair.h"
#include <cmath>
#include <cstring>
using namespace LAMMPS_NS;
#define BUFFACTOR 1.5
#define BUFMIN 1024
#define BIG 1.0e20
CommBrick::CommBrick(LAMMPS *lmp)
    : Comm(lmp), sendnum(nullptr), recvnum(nullptr), sendproc(nullptr),
      recvproc(nullptr), size_forward_recv(nullptr), size_reverse_send(nullptr),
      size_reverse_recv(nullptr), slablo(nullptr), slabhi(nullptr),
      multilo(nullptr), multihi(nullptr), multioldlo(nullptr),
      multioldhi(nullptr), cutghostmulti(nullptr), cutghostmultiold(nullptr),
      pbc_flag(nullptr), pbc(nullptr), firstrecv(nullptr), sendlist(nullptr),
      localsendlist(nullptr), maxsendlist(nullptr), buf_send(nullptr),
      buf_recv(nullptr) {
  style = Comm::BRICK;
  layout = Comm::LAYOUT_UNIFORM;
  pbc_flag = nullptr;
  init_buffers();
}
CommBrick::~CommBrick() {
  CommBrick::free_swap();
  if (sendlist)
    for (int i = 0; i < maxswap; i++)
      memory->destroy(sendlist[i]);
  if (localsendlist)
    memory->destroy(localsendlist);
  memory->sfree(sendlist);
  memory->destroy(maxsendlist);
  memory->destroy(buf_send);
  memory->destroy(buf_recv);
}
void CommBrick::init_buffers() {
  multilo = multihi = nullptr;
  cutghostmulti = nullptr;
  multioldlo = multioldhi = nullptr;
  cutghostmultiold = nullptr;
  buf_send = buf_recv = nullptr;
  maxsend = maxrecv = BUFMIN;
  CommBrick::grow_send(maxsend, 2);
  memory->create(buf_recv, maxrecv, "comm:buf_recv");
  nswap = 0;
  maxswap = 6;
  CommBrick::allocate_swap(maxswap);
  sendlist = (int **)memory->smalloc(maxswap * sizeof(int *), "comm:sendlist");
  memory->create(maxsendlist, maxswap, "comm:maxsendlist");
  for (int i = 0; i < maxswap; i++) {
    maxsendlist[i] = BUFMIN;
    memory->create(sendlist[i], BUFMIN, "comm:sendlist[i]");
  }
}
void CommBrick::init() {
  Comm::init();
  int bufextra_old = bufextra;
  init_exchange();
  if (bufextra > bufextra_old)
    grow_send(maxsend + bufextra, 2);
}
void CommBrick::setup() {
  int i, j;
  int ntypes = atom->ntypes;
  double *prd, *sublo, *subhi;
  double cut = get_comm_cutoff();
  if ((cut == 0.0) && (me == 0))
    error->warning(FLERR, "Communication cutoff is 0.0. No ghost atoms "
                          "will be generated. Atoms may get lost.");
  prd = domain->prd;
  sublo = domain->sublo;
  subhi = domain->subhi;
  cutghost[0] = cutghost[1] = cutghost[2] = cut;
  int *periodicity = domain->periodicity;
  int left, right;
  if (layout == Comm::LAYOUT_UNIFORM) {
    maxneed[0] = static_cast<int>(cutghost[0] * procgrid[0] / prd[0]) + 1;
    maxneed[1] = static_cast<int>(cutghost[1] * procgrid[1] / prd[1]) + 1;
    maxneed[2] = static_cast<int>(cutghost[2] * procgrid[2] / prd[2]) + 1;
    recvneed[0][0] = recvneed[0][1] = sendneed[0][0] = sendneed[0][1] =
        maxneed[0];
    recvneed[1][0] = recvneed[1][1] = sendneed[1][0] = sendneed[1][1] =
        maxneed[1];
    recvneed[2][0] = recvneed[2][1] = sendneed[2][0] = sendneed[2][1] =
        maxneed[2];
  }
  nswap = 2 * (maxneed[0] + maxneed[1] + maxneed[2]);
  if (nswap > maxswap)
    grow_swap(nswap);
  int dim, ineed;
  int iswap = 0;
  for (dim = 0; dim < 3; dim++) {
    for (ineed = 0; ineed < 2 * maxneed[dim]; ineed++) {
      pbc_flag[iswap] = 0;
      pbc[iswap][0] = pbc[iswap][1] = pbc[iswap][2] = pbc[iswap][3] =
          pbc[iswap][4] = pbc[iswap][5] = 0;
      if (ineed % 2 == 0) {
        sendproc[iswap] = procneigh[dim][0];
        recvproc[iswap] = procneigh[dim][1];
        if (mode == Comm::SINGLE) {
          if (ineed < 2)
            slablo[iswap] = -BIG;
          else
            slablo[iswap] = 0.5 * (sublo[dim] + subhi[dim]);
          slabhi[iswap] = sublo[dim] + cutghost[dim];
        }
        if (myloc[dim] == 0) {
          pbc_flag[iswap] = 1;
          pbc[iswap][dim] = 1;
          if (triclinic) {
            if (dim == 1)
              pbc[iswap][5] = 1;
            else if (dim == 2)
              pbc[iswap][4] = pbc[iswap][3] = 1;
          }
        }
      } else {
        sendproc[iswap] = procneigh[dim][1];
        recvproc[iswap] = procneigh[dim][0];
        if (mode == Comm::SINGLE) {
          slablo[iswap] = subhi[dim] - cutghost[dim];
          if (ineed < 2)
            slabhi[iswap] = BIG;
          else
            slabhi[iswap] = 0.5 * (sublo[dim] + subhi[dim]);
        } else {
          for (i = 1; i <= ntypes; i++) {
            multioldlo[iswap][i] = subhi[dim] - cutghostmultiold[i][dim];
            if (ineed < 2)
              multioldhi[iswap][i] = BIG;
            else
              multioldhi[iswap][i] = 0.5 * (sublo[dim] + subhi[dim]);
          }
        }
        if (myloc[dim] == procgrid[dim] - 1) {
          pbc_flag[iswap] = 1;
          pbc[iswap][dim] = -1;
          if (triclinic) {
            if (dim == 1)
              pbc[iswap][5] = -1;
            else if (dim == 2)
              pbc[iswap][4] = pbc[iswap][3] = -1;
          }
        }
      }
      iswap++;
    }
  }
}
void CommBrick::forward_comm(int) {
  int n;
  MPI_Request request;
  AtomVec *avec = atom->avec;
  double **x = atom->x;
  double *buf;
  for (int iswap = 0; iswap < nswap; iswap++) {
    if (sendproc[iswap] != me) {
      if (comm_x_only) {
        if (size_forward_recv[iswap]) {
          buf = x[firstrecv[iswap]];
          MPI_Irecv(buf, size_forward_recv[iswap], MPI_DOUBLE, recvproc[iswap],
                    0, world, &request);
        }
        n = avec->pack_comm(sendnum[iswap], sendlist[iswap], buf_send,
                            pbc_flag[iswap], pbc[iswap]);
        if (n)
          MPI_Send(buf_send, n, MPI_DOUBLE, sendproc[iswap], 0, world);
        if (size_forward_recv[iswap])
          MPI_Wait(&request, MPI_STATUS_IGNORE);
      } else if (ghost_velocity) {
        if (size_forward_recv[iswap])
          MPI_Irecv(buf_recv, size_forward_recv[iswap], MPI_DOUBLE,
                    recvproc[iswap], 0, world, &request);
        n = avec->pack_comm_vel(sendnum[iswap], sendlist[iswap], buf_send,
                                pbc_flag[iswap], pbc[iswap]);
        if (n)
          MPI_Send(buf_send, n, MPI_DOUBLE, sendproc[iswap], 0, world);
        if (size_forward_recv[iswap])
          MPI_Wait(&request, MPI_STATUS_IGNORE);
        avec->unpack_comm_vel(recvnum[iswap], firstrecv[iswap], buf_recv);
      } else {
        if (size_forward_recv[iswap])
          MPI_Irecv(buf_recv, size_forward_recv[iswap], MPI_DOUBLE,
                    recvproc[iswap], 0, world, &request);
        n = avec->pack_comm(sendnum[iswap], sendlist[iswap], buf_send,
                            pbc_flag[iswap], pbc[iswap]);
        if (n)
          MPI_Send(buf_send, n, MPI_DOUBLE, sendproc[iswap], 0, world);
        if (size_forward_recv[iswap])
          MPI_Wait(&request, MPI_STATUS_IGNORE);
        avec->unpack_comm(recvnum[iswap], firstrecv[iswap], buf_recv);
      }
    } else {
      if (comm_x_only) {
        if (sendnum[iswap])
          avec->pack_comm(sendnum[iswap], sendlist[iswap], x[firstrecv[iswap]],
                          pbc_flag[iswap], pbc[iswap]);
      } else if (ghost_velocity) {
        avec->pack_comm_vel(sendnum[iswap], sendlist[iswap], buf_send,
                            pbc_flag[iswap], pbc[iswap]);
        avec->unpack_comm_vel(recvnum[iswap], firstrecv[iswap], buf_send);
      } else {
        avec->pack_comm(sendnum[iswap], sendlist[iswap], buf_send,
                        pbc_flag[iswap], pbc[iswap]);
        avec->unpack_comm(recvnum[iswap], firstrecv[iswap], buf_send);
      }
    }
  }
}
void CommBrick::reverse_comm() {
  int n;
  MPI_Request request;
  AtomVec *avec = atom->avec;
  double **f = atom->f;
  double *buf;
  for (int iswap = nswap - 1; iswap >= 0; iswap--) {
    if (sendproc[iswap] != me) {
      if (comm_f_only) {
        if (size_reverse_recv[iswap])
          MPI_Irecv(buf_recv, size_reverse_recv[iswap], MPI_DOUBLE,
                    sendproc[iswap], 0, world, &request);
        if (size_reverse_send[iswap]) {
          buf = f[firstrecv[iswap]];
          MPI_Send(buf, size_reverse_send[iswap], MPI_DOUBLE, recvproc[iswap],
                   0, world);
        }
        if (size_reverse_recv[iswap])
          MPI_Wait(&request, MPI_STATUS_IGNORE);
      } else {
        if (size_reverse_recv[iswap])
          MPI_Irecv(buf_recv, size_reverse_recv[iswap], MPI_DOUBLE,
                    sendproc[iswap], 0, world, &request);
        n = avec->pack_reverse(recvnum[iswap], firstrecv[iswap], buf_send);
        if (n)
          MPI_Send(buf_send, n, MPI_DOUBLE, recvproc[iswap], 0, world);
        if (size_reverse_recv[iswap])
          MPI_Wait(&request, MPI_STATUS_IGNORE);
      }
      avec->unpack_reverse(sendnum[iswap], sendlist[iswap], buf_recv);
    } else {
      if (comm_f_only) {
        if (sendnum[iswap])
          avec->unpack_reverse(sendnum[iswap], sendlist[iswap],
                               f[firstrecv[iswap]]);
      } else {
        avec->pack_reverse(recvnum[iswap], firstrecv[iswap], buf_send);
        avec->unpack_reverse(sendnum[iswap], sendlist[iswap], buf_send);
      }
    }
  }
}
void CommBrick::exchange() {
  int i, m, nsend, nrecv, nrecv1, nrecv2, nlocal;
  double lo, hi, value;
  double **x;
  double *sublo, *subhi;
  MPI_Request request;
  AtomVec *avec = atom->avec;
  if (map_style != Atom::MAP_NONE)
    atom->map_clear();
  atom->nghost = 0;
  atom->avec->clear_bonus();
  if (maxexchange_fix_dynamic) {
    int bufextra_old = bufextra;
    init_exchange();
    if (bufextra > bufextra_old)
      grow_send(maxsend + bufextra, 2);
  }
  if (triclinic == 0) {
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }
  int dimension = domain->dimension;
  for (int dim = 0; dim < dimension; dim++) {
    x = atom->x;
    lo = sublo[dim];
    hi = subhi[dim];
    nlocal = atom->nlocal;
    i = nsend = 0;
    while (i < nlocal) {
      if (x[i][dim] < lo || x[i][dim] >= hi) {
        if (nsend > maxsend)
          grow_send(nsend, 1);
        nsend += avec->pack_exchange(i, &buf_send[nsend]);
        avec->copy(nlocal - 1, i, 1);
        nlocal--;
      } else
        i++;
    }
    atom->nlocal = nlocal;
    if (procgrid[dim] == 1)
      nrecv = 0;
    else {
      MPI_Sendrecv(&nsend, 1, MPI_INT, procneigh[dim][0], 0, &nrecv1, 1,
                   MPI_INT, procneigh[dim][1], 0, world, MPI_STATUS_IGNORE);
      nrecv = nrecv1;
      if (procgrid[dim] > 2) {
        MPI_Sendrecv(&nsend, 1, MPI_INT, procneigh[dim][1], 0, &nrecv2, 1,
                     MPI_INT, procneigh[dim][0], 0, world, MPI_STATUS_IGNORE);
        nrecv += nrecv2;
      }
      if (nrecv > maxrecv)
        grow_recv(nrecv);
      MPI_Irecv(buf_recv, nrecv1, MPI_DOUBLE, procneigh[dim][1], 0, world,
                &request);
      MPI_Send(buf_send, nsend, MPI_DOUBLE, procneigh[dim][0], 0, world);
      MPI_Wait(&request, MPI_STATUS_IGNORE);
      if (procgrid[dim] > 2) {
        MPI_Irecv(&buf_recv[nrecv1], nrecv2, MPI_DOUBLE, procneigh[dim][0], 0,
                  world, &request);
        MPI_Send(buf_send, nsend, MPI_DOUBLE, procneigh[dim][1], 0, world);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
      }
    }
    m = 0;
    while (m < nrecv) {
      value = buf_recv[m + dim + 1];
      if (value >= lo && value < hi)
        m += avec->unpack_exchange(&buf_recv[m]);
      else
        m += static_cast<int>(buf_recv[m]);
    }
  }
  if (atom->firstgroupname)
    atom->first_reorder();
}
void CommBrick::borders() {
  int i, n, itype, icollection, iswap, dim, ineed, twoneed;
  int nsend, nrecv, sendflag, nfirst, nlast, ngroup, nprior;
  double lo, hi;
  int *type;
  int *collection;
  double **x;
  double *buf, *mlo, *mhi;
  MPI_Request request;
  AtomVec *avec = atom->avec;
  iswap = 0;
  smax = rmax = 0;
  for (dim = 0; dim < 3; dim++) {
    nlast = 0;
    twoneed = 2 * maxneed[dim];
    for (ineed = 0; ineed < twoneed; ineed++) {
      x = atom->x;
      if (mode == Comm::SINGLE) {
        lo = slablo[iswap];
        hi = slabhi[iswap];
      } else {
        type = atom->type;
        mlo = multioldlo[iswap];
        mhi = multioldhi[iswap];
      }
      if (ineed % 2 == 0) {
        nfirst = nlast;
        nlast = atom->nlocal + atom->nghost;
      }
      nsend = 0;
      if (ineed / 2 >= sendneed[dim][ineed % 2])
        sendflag = 0;
      else
        sendflag = 1;
      if (sendflag) {
        if (!bordergroup || ineed >= 2) {
          if (mode == Comm::SINGLE) {
            for (i = nfirst; i < nlast; i++)
              if (x[i][dim] >= lo && x[i][dim] <= hi) {
                if (nsend == maxsendlist[iswap])
                  grow_list(iswap, nsend);
                sendlist[iswap][nsend++] = i;
              }
          } else {
            for (i = nfirst; i < nlast; i++) {
              itype = type[i];
              if (x[i][dim] >= mlo[itype] && x[i][dim] <= mhi[itype]) {
                if (nsend == maxsendlist[iswap])
                  grow_list(iswap, nsend);
                sendlist[iswap][nsend++] = i;
              }
            }
          }
        } else {
          if (mode == Comm::SINGLE) {
            ngroup = atom->nfirst;
            for (i = 0; i < ngroup; i++)
              if (x[i][dim] >= lo && x[i][dim] <= hi) {
                if (nsend == maxsendlist[iswap])
                  grow_list(iswap, nsend);
                sendlist[iswap][nsend++] = i;
              }
            for (i = atom->nlocal; i < nlast; i++)
              if (x[i][dim] >= lo && x[i][dim] <= hi) {
                if (nsend == maxsendlist[iswap])
                  grow_list(iswap, nsend);
                sendlist[iswap][nsend++] = i;
              }
          } else {
            ngroup = atom->nfirst;
            for (i = 0; i < ngroup; i++) {
              itype = type[i];
              if (x[i][dim] >= mlo[itype] && x[i][dim] <= mhi[itype]) {
                if (nsend == maxsendlist[iswap])
                  grow_list(iswap, nsend);
                sendlist[iswap][nsend++] = i;
              }
            }
            for (i = atom->nlocal; i < nlast; i++) {
              itype = type[i];
              if (x[i][dim] >= mlo[itype] && x[i][dim] <= mhi[itype]) {
                if (nsend == maxsendlist[iswap])
                  grow_list(iswap, nsend);
                sendlist[iswap][nsend++] = i;
              }
            }
          }
        }
      }
      if (nsend * size_border > maxsend)
        grow_send(nsend * size_border, 0);
      if (ghost_velocity)
        n = avec->pack_border_vel(nsend, sendlist[iswap], buf_send,
                                  pbc_flag[iswap], pbc[iswap]);
      else
        n = avec->pack_border(nsend, sendlist[iswap], buf_send, pbc_flag[iswap],
                              pbc[iswap]);
      if (sendproc[iswap] != me) {
        MPI_Sendrecv(&nsend, 1, MPI_INT, sendproc[iswap], 0, &nrecv, 1, MPI_INT,
                     recvproc[iswap], 0, world, MPI_STATUS_IGNORE);
        if (nrecv * size_border > maxrecv)
          grow_recv(nrecv * size_border);
        if (nrecv)
          MPI_Irecv(buf_recv, nrecv * size_border, MPI_DOUBLE, recvproc[iswap],
                    0, world, &request);
        if (n)
          MPI_Send(buf_send, n, MPI_DOUBLE, sendproc[iswap], 0, world);
        if (nrecv)
          MPI_Wait(&request, MPI_STATUS_IGNORE);
        buf = buf_recv;
      } else {
        nrecv = nsend;
        buf = buf_send;
      }
      if (ghost_velocity)
        avec->unpack_border_vel(nrecv, atom->nlocal + atom->nghost, buf);
      else
        avec->unpack_border(nrecv, atom->nlocal + atom->nghost, buf);
      smax = MAX(smax, nsend);
      rmax = MAX(rmax, nrecv);
      sendnum[iswap] = nsend;
      recvnum[iswap] = nrecv;
      size_forward_recv[iswap] = nrecv * size_forward;
      size_reverse_send[iswap] = nrecv * size_reverse;
      size_reverse_recv[iswap] = nsend * size_reverse;
      firstrecv[iswap] = atom->nlocal + atom->nghost;
      nprior = atom->nlocal + atom->nghost;
      atom->nghost += nrecv;
      iswap++;
    }
  }
  int max = MAX(maxforward * smax, maxreverse * rmax);
  if (max > maxsend)
    grow_send(max, 0);
  max = MAX(maxforward * rmax, maxreverse * smax);
  if (max > maxrecv)
    grow_recv(max);
  if (map_style != Atom::MAP_NONE)
    atom->map_set();
}
void CommBrick::forward_comm(Pair *pair) {
  int iswap, n;
  double *buf;
  MPI_Request request;
  int nsize = pair->comm_forward;
  for (iswap = 0; iswap < nswap; iswap++) {
    n = pair->pack_forward_comm(sendnum[iswap], sendlist[iswap], buf_send,
                                pbc_flag[iswap], pbc[iswap]);
    if (sendproc[iswap] != me) {
      if (recvnum[iswap])
        MPI_Irecv(buf_recv, nsize * recvnum[iswap], MPI_DOUBLE, recvproc[iswap],
                  0, world, &request);
      if (sendnum[iswap])
        MPI_Send(buf_send, n, MPI_DOUBLE, sendproc[iswap], 0, world);
      if (recvnum[iswap])
        MPI_Wait(&request, MPI_STATUS_IGNORE);
      buf = buf_recv;
    } else
      buf = buf_send;
    pair->unpack_forward_comm(recvnum[iswap], firstrecv[iswap], buf);
  }
}
void CommBrick::reverse_comm(Pair *pair) {
  int iswap, n;
  double *buf;
  MPI_Request request;
  int nsize = MAX(pair->comm_reverse, pair->comm_reverse_off);
  for (iswap = nswap - 1; iswap >= 0; iswap--) {
    n = pair->pack_reverse_comm(recvnum[iswap], firstrecv[iswap], buf_send);
    if (sendproc[iswap] != me) {
      if (sendnum[iswap])
        MPI_Irecv(buf_recv, nsize * sendnum[iswap], MPI_DOUBLE, sendproc[iswap],
                  0, world, &request);
      if (recvnum[iswap])
        MPI_Send(buf_send, n, MPI_DOUBLE, recvproc[iswap], 0, world);
      if (sendnum[iswap])
        MPI_Wait(&request, MPI_STATUS_IGNORE);
      buf = buf_recv;
    } else
      buf = buf_send;
    pair->unpack_reverse_comm(sendnum[iswap], sendlist[iswap], buf);
  }
}
void CommBrick::forward_comm(Fix *fix, int size) {
  int iswap, n, nsize;
  double *buf;
  MPI_Request request;
  if (size)
    nsize = size;
  else
    nsize = fix->comm_forward;
  for (iswap = 0; iswap < nswap; iswap++) {
    n = fix->pack_forward_comm(sendnum[iswap], sendlist[iswap], buf_send,
                               pbc_flag[iswap], pbc[iswap]);
    if (sendproc[iswap] != me) {
      if (recvnum[iswap])
        MPI_Irecv(buf_recv, nsize * recvnum[iswap], MPI_DOUBLE, recvproc[iswap],
                  0, world, &request);
      if (sendnum[iswap])
        MPI_Send(buf_send, n, MPI_DOUBLE, sendproc[iswap], 0, world);
      if (recvnum[iswap])
        MPI_Wait(&request, MPI_STATUS_IGNORE);
      buf = buf_recv;
    } else
      buf = buf_send;
    fix->unpack_forward_comm(recvnum[iswap], firstrecv[iswap], buf);
  }
}
void CommBrick::reverse_comm(Fix *fix, int size) {
  int iswap, n, nsize;
  double *buf;
  MPI_Request request;
  if (size)
    nsize = size;
  else
    nsize = fix->comm_reverse;
  for (iswap = nswap - 1; iswap >= 0; iswap--) {
    n = fix->pack_reverse_comm(recvnum[iswap], firstrecv[iswap], buf_send);
    if (sendproc[iswap] != me) {
      if (sendnum[iswap])
        MPI_Irecv(buf_recv, nsize * sendnum[iswap], MPI_DOUBLE, sendproc[iswap],
                  0, world, &request);
      if (recvnum[iswap])
        MPI_Send(buf_send, n, MPI_DOUBLE, recvproc[iswap], 0, world);
      if (sendnum[iswap])
        MPI_Wait(&request, MPI_STATUS_IGNORE);
      buf = buf_recv;
    } else
      buf = buf_send;
    fix->unpack_reverse_comm(sendnum[iswap], sendlist[iswap], buf);
  }
}
void CommBrick::reverse_comm_variable(Fix *fix) {
  int iswap, nsend, nrecv;
  double *buf;
  MPI_Request request;
  for (iswap = nswap - 1; iswap >= 0; iswap--) {
    nsend = fix->pack_reverse_comm_size(recvnum[iswap], firstrecv[iswap]);
    if (nsend > maxsend)
      grow_send(nsend, 0);
    nsend = fix->pack_reverse_comm(recvnum[iswap], firstrecv[iswap], buf_send);
    if (sendproc[iswap] != me) {
      MPI_Sendrecv(&nsend, 1, MPI_INT, recvproc[iswap], 0, &nrecv, 1, MPI_INT,
                   sendproc[iswap], 0, world, MPI_STATUS_IGNORE);
      if (sendnum[iswap]) {
        if (nrecv > maxrecv)
          grow_recv(nrecv);
        MPI_Irecv(buf_recv, maxrecv, MPI_DOUBLE, sendproc[iswap], 0, world,
                  &request);
      }
      if (recvnum[iswap])
        MPI_Send(buf_send, nsend, MPI_DOUBLE, recvproc[iswap], 0, world);
      if (sendnum[iswap])
        MPI_Wait(&request, MPI_STATUS_IGNORE);
      buf = buf_recv;
    } else
      buf = buf_send;
    fix->unpack_reverse_comm(sendnum[iswap], sendlist[iswap], buf);
  }
}
void CommBrick::grow_send(int n, int flag) {
  if (flag == 0) {
    maxsend = static_cast<int>(BUFFACTOR * n);
    memory->destroy(buf_send);
    memory->create(buf_send, maxsend + bufextra, "comm:buf_send");
  } else if (flag == 1) {
    maxsend = static_cast<int>(BUFFACTOR * n);
    memory->grow(buf_send, maxsend + bufextra, "comm:buf_send");
  } else {
    memory->destroy(buf_send);
    memory->grow(buf_send, maxsend + bufextra, "comm:buf_send");
  }
}
void CommBrick::grow_recv(int n) {
  maxrecv = static_cast<int>(BUFFACTOR * n);
  memory->destroy(buf_recv);
  memory->create(buf_recv, maxrecv, "comm:buf_recv");
}
void CommBrick::grow_list(int iswap, int n) {
  maxsendlist[iswap] = static_cast<int>(BUFFACTOR * n);
  memory->grow(sendlist[iswap], maxsendlist[iswap], "comm:sendlist[iswap]");
}
void CommBrick::grow_swap(int n) {
  free_swap();
  allocate_swap(n);
  sendlist =
      (int **)memory->srealloc(sendlist, n * sizeof(int *), "comm:sendlist");
  memory->grow(maxsendlist, n, "comm:maxsendlist");
  for (int i = maxswap; i < n; i++) {
    maxsendlist[i] = BUFMIN;
    memory->create(sendlist[i], BUFMIN, "comm:sendlist[i]");
  }
  maxswap = n;
}
void CommBrick::allocate_swap(int n) {
  memory->create(sendnum, n, "comm:sendnum");
  memory->create(recvnum, n, "comm:recvnum");
  memory->create(sendproc, n, "comm:sendproc");
  memory->create(recvproc, n, "comm:recvproc");
  memory->create(size_forward_recv, n, "comm:size");
  memory->create(size_reverse_send, n, "comm:size");
  memory->create(size_reverse_recv, n, "comm:size");
  memory->create(slablo, n, "comm:slablo");
  memory->create(slabhi, n, "comm:slabhi");
  memory->create(firstrecv, n, "comm:firstrecv");
  memory->create(pbc_flag, n, "comm:pbc_flag");
  memory->create(pbc, n, 6, "comm:pbc");
}
void CommBrick::free_swap() {
  memory->destroy(sendnum);
  memory->destroy(recvnum);
  memory->destroy(sendproc);
  memory->destroy(recvproc);
  memory->destroy(size_forward_recv);
  memory->destroy(size_reverse_send);
  memory->destroy(size_reverse_recv);
  memory->destroy(slablo);
  memory->destroy(slabhi);
  memory->destroy(firstrecv);
  memory->destroy(pbc_flag);
  memory->destroy(pbc);
}
