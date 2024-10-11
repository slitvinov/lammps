#include "comm_tiled.h"
#include "atom.h"
#include "atom_vec.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "memory.h"
#include "neighbor.h"
#include "pair.h"
#include <cmath>
#include <cstring>
#include <float.h>
using namespace LAMMPS_NS;
#define BUFFACTOR 1.5
#define BUFFACTOR 1.5
#define BUFMIN 1024
#define EPSILON 1.0e-6
#define DELTA_PROCS 16
CommTiled::CommTiled(LAMMPS *lmp) : Comm(lmp)
{
  style = Comm::TILED;
  layout = Comm::LAYOUT_UNIFORM;
  pbc_flag = nullptr;
  buf_send = nullptr;
  buf_recv = nullptr;
  overlap = nullptr;
  rcbinfo = nullptr;
  cutghostmulti = nullptr;
  cutghostmultiold = nullptr;
  init_buffers();
}
CommTiled::CommTiled(LAMMPS * , Comm *oldcomm) : Comm(*oldcomm)
{
  style = Comm::TILED;
  layout = oldcomm->layout;
  Comm::copy_arrays(oldcomm);
  init_buffers();
}
CommTiled::~CommTiled()
{
  memory->destroy(buf_send);
  memory->destroy(buf_recv);
  memory->destroy(overlap);
  deallocate_swap(maxswap);
  memory->sfree(rcbinfo);
  memory->destroy(cutghostmulti);
  memory->destroy(cutghostmultiold);
}
void CommTiled::init_buffers()
{
  buf_send = buf_recv = nullptr;
  maxsend = maxrecv = BUFMIN;
  grow_send(maxsend, 2);
  memory->create(buf_recv, maxrecv, "comm:buf_recv");
  maxoverlap = 0;
  overlap = nullptr;
  rcbinfo = nullptr;
  cutghostmulti = nullptr;
  cutghostmultiold = nullptr;
  sendbox_multi = nullptr;
  sendbox_multiold = nullptr;
  maxswap = 6;
  allocate_swap(maxswap);
}
void CommTiled::init()
{
  Comm::init();
  nswap = 2 * domain->dimension;
  memory->destroy(cutghostmulti);
  if (mode == Comm::MULTI) {
    if (ncollections != neighbor->ncollections) { ncollections = neighbor->ncollections; }
    if (cutusermulti && ncollections != ncollections_cutoff) {
      if (me == 0)
        error->warning(FLERR,
                       "cutoff/multi settings discarded, must be defined"
                       " after customizing collections in neigh_modify");
      memory->destroy(cutusermulti);
      cutusermulti = nullptr;
    }
    for (int i = 0; i < maxswap; i++) grow_swap_send_multi(i, DELTA_PROCS);
    memory->create(cutghostmulti, ncollections, 3, "comm:cutghostmulti");
  }
  memory->destroy(cutghostmultiold);
  if (mode == Comm::MULTIOLD)
    memory->create(cutghostmultiold, atom->ntypes + 1, 3, "comm:cutghostmultiold");
  int bufextra_old = bufextra;
  init_exchange();
  if (bufextra > bufextra_old) grow_send(maxsend + bufextra, 2);
}
void CommTiled::setup()
{
  int i, j, n;
  dimension = domain->dimension;
  int *periodicity = domain->periodicity;
  int ntypes = atom->ntypes;
  if (triclinic == 0) {
    prd = domain->prd;
    boxlo = domain->boxlo;
    boxhi = domain->boxhi;
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    prd = domain->prd_lamda;
    boxlo = domain->boxlo_lamda;
    boxhi = domain->boxhi_lamda;
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }
  if (layout != Comm::LAYOUT_TILED) {
    box_drop = &CommTiled::box_drop_brick;
    box_other = &CommTiled::box_other_brick;
    box_touch = &CommTiled::box_touch_brick;
    point_drop = &CommTiled::point_drop_brick;
  } else {
    box_drop = &CommTiled::box_drop_tiled;
    box_other = &CommTiled::box_other_tiled;
    box_touch = &CommTiled::box_touch_tiled;
    point_drop = &CommTiled::point_drop_tiled;
  }
  if (layout == Comm::LAYOUT_TILED) coord2proc_setup();
  if (mode == Comm::MULTI) {
    double **cutcollectionsq = neighbor->cutcollectionsq;
    neighbor->build_collection(0);
    for (i = 0; i < ncollections; i++) {
      if (cutusermulti) {
        cutghostmulti[i][0] = cutusermulti[i];
        cutghostmulti[i][1] = cutusermulti[i];
        cutghostmulti[i][2] = cutusermulti[i];
      } else {
        cutghostmulti[i][0] = 0.0;
        cutghostmulti[i][1] = 0.0;
        cutghostmulti[i][2] = 0.0;
      }
      for (j = 0; j < ncollections; j++) {
        if (multi_reduce && (cutcollectionsq[j][j] > cutcollectionsq[i][i])) continue;
        cutghostmulti[i][0] = MAX(cutghostmulti[i][0], sqrt(cutcollectionsq[i][j]));
        cutghostmulti[i][1] = MAX(cutghostmulti[i][1], sqrt(cutcollectionsq[i][j]));
        cutghostmulti[i][2] = MAX(cutghostmulti[i][2], sqrt(cutcollectionsq[i][j]));
      }
    }
  }
  if (mode == Comm::MULTIOLD) {
    double *cuttype = neighbor->cuttype;
    for (i = 1; i <= ntypes; i++) {
      double tmp = 0.0;
      if (cutusermultiold) tmp = cutusermultiold[i];
      cutghostmultiold[i][0] = MAX(tmp, cuttype[i]);
      cutghostmultiold[i][1] = MAX(tmp, cuttype[i]);
      cutghostmultiold[i][2] = MAX(tmp, cuttype[i]);
    }
  }
  double cut = get_comm_cutoff();
  if ((cut == 0.0) && (me == 0))
    error->warning(FLERR,
                   "Communication cutoff is 0.0. No ghost atoms "
                   "will be generated. Atoms may get lost.");
  if (triclinic == 0)
    cutghost[0] = cutghost[1] = cutghost[2] = cut;
  else {
    double *h_inv = domain->h_inv;
    double length0, length1, length2;
    length0 = sqrt(h_inv[0] * h_inv[0] + h_inv[5] * h_inv[5] + h_inv[4] * h_inv[4]);
    cutghost[0] = cut * length0;
    length1 = sqrt(h_inv[1] * h_inv[1] + h_inv[3] * h_inv[3]);
    cutghost[1] = cut * length1;
    length2 = h_inv[2];
    cutghost[2] = cut * length2;
    if (mode == Comm::MULTI) {
      for (i = 0; i < ncollections; i++) {
        cutghostmulti[i][0] *= length0;
        cutghostmulti[i][1] *= length1;
        cutghostmulti[i][2] *= length2;
      }
    }
    if (mode == Comm::MULTIOLD) {
      for (i = 1; i <= ntypes; i++) {
        cutghostmultiold[i][0] *= length0;
        cutghostmultiold[i][1] *= length1;
        cutghostmultiold[i][2] *= length2;
      }
    }
  }
  if ((periodicity[0] && cutghost[0] > prd[0]) || (periodicity[1] && cutghost[1] > prd[1]) ||
      (dimension == 3 && periodicity[2] && cutghost[2] > prd[2]))
    error->all(FLERR,
               "Communication cutoff for comm_style tiled "
               "cannot exceed periodic box length");
  int cutzero = 0;
  if (cut == 0.0) {
    cutzero = 1;
    cut = MIN(prd[0], prd[1]);
    if (dimension == 3) cut = MIN(cut, prd[2]);
    cut *= EPSILON * EPSILON;
    cutghost[0] = cutghost[1] = cutghost[2] = cut;
  }
  int noverlap1, indexme;
  double lo1[3], hi1[3], lo2[3], hi2[3];
  int one, two;
  int iswap = 0;
  for (int idim = 0; idim < dimension; idim++) {
    for (int idir = 0; idir < 2; idir++) {
      one = 1;
      lo1[0] = sublo[0];
      lo1[1] = sublo[1];
      lo1[2] = sublo[2];
      hi1[0] = subhi[0];
      hi1[1] = subhi[1];
      hi1[2] = subhi[2];
      if (idir == 0) {
        lo1[idim] = sublo[idim] - cutghost[idim];
        hi1[idim] = sublo[idim];
      } else {
        lo1[idim] = subhi[idim];
        hi1[idim] = subhi[idim] + cutghost[idim];
      }
      two = 0;
      if (idir == 0 && periodicity[idim] && lo1[idim] < boxlo[idim]) two = 1;
      if (idir == 1 && periodicity[idim] && hi1[idim] > boxhi[idim]) two = 1;
      if (two) {
        lo2[0] = sublo[0];
        lo2[1] = sublo[1];
        lo2[2] = sublo[2];
        hi2[0] = subhi[0];
        hi2[1] = subhi[1];
        hi2[2] = subhi[2];
        if (idir == 0) {
          lo2[idim] = lo1[idim] + prd[idim];
          hi2[idim] = boxhi[idim];
          if (sublo[idim] == boxlo[idim]) one = 0;
        } else {
          lo2[idim] = boxlo[idim];
          hi2[idim] = hi1[idim] - prd[idim];
          if (subhi[idim] == boxhi[idim]) one = 0;
        }
      }
      if (one) {
        if (idir == 0)
          lo1[idim] = MAX(lo1[idim], boxlo[idim]);
        else
          hi1[idim] = MIN(hi1[idim], boxhi[idim]);
        if (lo1[idim] == hi1[idim]) one = 0;
      }
      indexme = -1;
      noverlap = 0;
      if (one) (this->*box_drop)(idim, lo1, hi1, indexme);
      noverlap1 = noverlap;
      if (two) (this->*box_drop)(idim, lo2, hi2, indexme);
      if (indexme >= 0) {
        int tmp = overlap[noverlap - 1];
        overlap[noverlap - 1] = overlap[indexme];
        overlap[indexme] = tmp;
      }
      if (noverlap > nprocmax[iswap]) {
        int oldmax = nprocmax[iswap];
        while (nprocmax[iswap] < noverlap) nprocmax[iswap] += DELTA_PROCS;
        grow_swap_send(iswap, nprocmax[iswap], oldmax);
        if (idir == 0)
          grow_swap_recv(iswap + 1, nprocmax[iswap]);
        else
          grow_swap_recv(iswap - 1, nprocmax[iswap]);
      }
      if (noverlap && overlap[noverlap - 1] == me)
        sendself[iswap] = 1;
      else
        sendself[iswap] = 0;
      if (noverlap && noverlap - sendself[iswap])
        sendother[iswap] = 1;
      else
        sendother[iswap] = 0;
      nsendproc[iswap] = noverlap;
      for (i = 0; i < noverlap; i++) sendproc[iswap][i] = overlap[i];
      if (idir == 0) {
        recvother[iswap + 1] = sendother[iswap];
        nrecvproc[iswap + 1] = noverlap;
        for (i = 0; i < noverlap; i++) recvproc[iswap + 1][i] = overlap[i];
      } else {
        recvother[iswap - 1] = sendother[iswap];
        nrecvproc[iswap - 1] = noverlap;
        for (i = 0; i < noverlap; i++) recvproc[iswap - 1][i] = overlap[i];
      }
      double oboxlo[3], oboxhi[3], sbox[6], sbox_multi[6], sbox_multiold[6];
      if (mode == Comm::SINGLE) {
        for (i = 0; i < noverlap; i++) {
          pbc_flag[iswap][i] = 0;
          pbc[iswap][i][0] = pbc[iswap][i][1] = pbc[iswap][i][2] = pbc[iswap][i][3] =
              pbc[iswap][i][4] = pbc[iswap][i][5] = 0;
          (this->*box_other)(idim, idir, overlap[i], oboxlo, oboxhi);
          if (i < noverlap1) {
            sbox[0] = MAX(oboxlo[0], lo1[0]);
            sbox[1] = MAX(oboxlo[1], lo1[1]);
            sbox[2] = MAX(oboxlo[2], lo1[2]);
            sbox[3] = MIN(oboxhi[0], hi1[0]);
            sbox[4] = MIN(oboxhi[1], hi1[1]);
            sbox[5] = MIN(oboxhi[2], hi1[2]);
          } else {
            pbc_flag[iswap][i] = 1;
            if (idir == 0)
              pbc[iswap][i][idim] = 1;
            else
              pbc[iswap][i][idim] = -1;
            if (triclinic) {
              if (idim == 1) pbc[iswap][i][5] = pbc[iswap][i][idim];
              if (idim == 2) pbc[iswap][i][4] = pbc[iswap][i][3] = pbc[iswap][i][idim];
            }
            sbox[0] = MAX(oboxlo[0], lo2[0]);
            sbox[1] = MAX(oboxlo[1], lo2[1]);
            sbox[2] = MAX(oboxlo[2], lo2[2]);
            sbox[3] = MIN(oboxhi[0], hi2[0]);
            sbox[4] = MIN(oboxhi[1], hi2[1]);
            sbox[5] = MIN(oboxhi[2], hi2[2]);
          }
          if (idir == 0) {
            sbox[idim] = sublo[idim];
            if (i < noverlap1)
              sbox[3 + idim] = MIN(sbox[3 + idim] + cutghost[idim], subhi[idim]);
            else
              sbox[3 + idim] = MIN(sbox[3 + idim] - prd[idim] + cutghost[idim], subhi[idim]);
          } else {
            if (i < noverlap1)
              sbox[idim] = MAX(sbox[idim] - cutghost[idim], sublo[idim]);
            else
              sbox[idim] = MAX(sbox[idim] + prd[idim] - cutghost[idim], sublo[idim]);
            sbox[3 + idim] = subhi[idim];
          }
          if (idim >= 1) {
            if (sbox[0] == oboxlo[0]) sbox[0] -= cutghost[0];
            if (sbox[3] == oboxhi[0]) sbox[3] += cutghost[0];
          }
          if (idim == 2) {
            if (sbox[1] == oboxlo[1]) sbox[1] -= cutghost[1];
            if (sbox[4] == oboxhi[1]) sbox[4] += cutghost[1];
          }
          memcpy(sendbox[iswap][i], sbox, 6 * sizeof(double));
        }
      }
      if (mode == Comm::MULTI) {
        for (i = 0; i < noverlap; i++) {
          pbc_flag[iswap][i] = 0;
          pbc[iswap][i][0] = pbc[iswap][i][1] = pbc[iswap][i][2] = pbc[iswap][i][3] =
              pbc[iswap][i][4] = pbc[iswap][i][5] = 0;
          (this->*box_other)(idim, idir, overlap[i], oboxlo, oboxhi);
          if (i < noverlap1) {
            sbox[0] = MAX(oboxlo[0], lo1[0]);
            sbox[1] = MAX(oboxlo[1], lo1[1]);
            sbox[2] = MAX(oboxlo[2], lo1[2]);
            sbox[3] = MIN(oboxhi[0], hi1[0]);
            sbox[4] = MIN(oboxhi[1], hi1[1]);
            sbox[5] = MIN(oboxhi[2], hi1[2]);
          } else {
            pbc_flag[iswap][i] = 1;
            if (idir == 0)
              pbc[iswap][i][idim] = 1;
            else
              pbc[iswap][i][idim] = -1;
            if (triclinic) {
              if (idim == 1) pbc[iswap][i][5] = pbc[iswap][i][idim];
              if (idim == 2) pbc[iswap][i][4] = pbc[iswap][i][3] = pbc[iswap][i][idim];
            }
            sbox[0] = MAX(oboxlo[0], lo2[0]);
            sbox[1] = MAX(oboxlo[1], lo2[1]);
            sbox[2] = MAX(oboxlo[2], lo2[2]);
            sbox[3] = MIN(oboxhi[0], hi2[0]);
            sbox[4] = MIN(oboxhi[1], hi2[1]);
            sbox[5] = MIN(oboxhi[2], hi2[2]);
          }
          for (int icollection = 0; icollection < ncollections; icollection++) {
            sbox_multi[0] = sbox[0];
            sbox_multi[1] = sbox[1];
            sbox_multi[2] = sbox[2];
            sbox_multi[3] = sbox[3];
            sbox_multi[4] = sbox[4];
            sbox_multi[5] = sbox[5];
            if (idir == 0) {
              sbox_multi[idim] = sublo[idim];
              if (i < noverlap1)
                sbox_multi[3 + idim] =
                    MIN(sbox_multi[3 + idim] + cutghostmulti[icollection][idim], subhi[idim]);
              else
                sbox_multi[3 + idim] =
                    MIN(sbox_multi[3 + idim] - prd[idim] + cutghostmulti[icollection][idim],
                        subhi[idim]);
            } else {
              if (i < noverlap1)
                sbox_multi[idim] =
                    MAX(sbox_multi[idim] - cutghostmulti[icollection][idim], sublo[idim]);
              else
                sbox_multi[idim] = MAX(
                    sbox_multi[idim] + prd[idim] - cutghostmulti[icollection][idim], sublo[idim]);
              sbox_multi[3 + idim] = subhi[idim];
            }
            if (idim >= 1) {
              if (sbox_multi[0] == oboxlo[0]) sbox_multi[0] -= cutghostmulti[icollection][idim];
              if (sbox_multi[3] == oboxhi[0]) sbox_multi[3] += cutghostmulti[icollection][idim];
            }
            if (idim == 2) {
              if (sbox_multi[1] == oboxlo[1]) sbox_multi[1] -= cutghostmulti[icollection][idim];
              if (sbox_multi[4] == oboxhi[1]) sbox_multi[4] += cutghostmulti[icollection][idim];
            }
            memcpy(sendbox_multi[iswap][i][icollection], sbox_multi, 6 * sizeof(double));
          }
        }
      }
      if (mode == Comm::MULTIOLD) {
        for (i = 0; i < noverlap; i++) {
          pbc_flag[iswap][i] = 0;
          pbc[iswap][i][0] = pbc[iswap][i][1] = pbc[iswap][i][2] = pbc[iswap][i][3] =
              pbc[iswap][i][4] = pbc[iswap][i][5] = 0;
          (this->*box_other)(idim, idir, overlap[i], oboxlo, oboxhi);
          if (i < noverlap1) {
            sbox[0] = MAX(oboxlo[0], lo1[0]);
            sbox[1] = MAX(oboxlo[1], lo1[1]);
            sbox[2] = MAX(oboxlo[2], lo1[2]);
            sbox[3] = MIN(oboxhi[0], hi1[0]);
            sbox[4] = MIN(oboxhi[1], hi1[1]);
            sbox[5] = MIN(oboxhi[2], hi1[2]);
          } else {
            pbc_flag[iswap][i] = 1;
            if (idir == 0)
              pbc[iswap][i][idim] = 1;
            else
              pbc[iswap][i][idim] = -1;
            if (triclinic) {
              if (idim == 1) pbc[iswap][i][5] = pbc[iswap][i][idim];
              if (idim == 2) pbc[iswap][i][4] = pbc[iswap][i][3] = pbc[iswap][i][idim];
            }
            sbox[0] = MAX(oboxlo[0], lo2[0]);
            sbox[1] = MAX(oboxlo[1], lo2[1]);
            sbox[2] = MAX(oboxlo[2], lo2[2]);
            sbox[3] = MIN(oboxhi[0], hi2[0]);
            sbox[4] = MIN(oboxhi[1], hi2[1]);
            sbox[5] = MIN(oboxhi[2], hi2[2]);
          }
          for (int itype = 1; itype <= atom->ntypes; itype++) {
            sbox_multiold[0] = sbox[0];
            sbox_multiold[1] = sbox[1];
            sbox_multiold[2] = sbox[2];
            sbox_multiold[3] = sbox[3];
            sbox_multiold[4] = sbox[4];
            sbox_multiold[5] = sbox[5];
            if (idir == 0) {
              sbox_multiold[idim] = sublo[idim];
              if (i < noverlap1)
                sbox_multiold[3 + idim] =
                    MIN(sbox_multiold[3 + idim] + cutghostmultiold[itype][idim], subhi[idim]);
              else
                sbox_multiold[3 + idim] =
                    MIN(sbox_multiold[3 + idim] - prd[idim] + cutghostmultiold[itype][idim],
                        subhi[idim]);
            } else {
              if (i < noverlap1)
                sbox_multiold[idim] =
                    MAX(sbox_multiold[idim] - cutghostmultiold[itype][idim], sublo[idim]);
              else
                sbox_multiold[idim] = MAX(
                    sbox_multiold[idim] + prd[idim] - cutghostmultiold[itype][idim], sublo[idim]);
              sbox_multiold[3 + idim] = subhi[idim];
            }
            if (idim >= 1) {
              if (sbox_multiold[0] == oboxlo[0]) sbox_multiold[0] -= cutghostmultiold[itype][idim];
              if (sbox_multiold[3] == oboxhi[0]) sbox_multiold[3] += cutghostmultiold[itype][idim];
            }
            if (idim == 2) {
              if (sbox_multiold[1] == oboxlo[1]) sbox_multiold[1] -= cutghostmultiold[itype][idim];
              if (sbox_multiold[4] == oboxhi[1]) sbox_multiold[4] += cutghostmultiold[itype][idim];
            }
            memcpy(sendbox_multiold[iswap][i][itype], sbox_multiold, 6 * sizeof(double));
          }
        }
      }
      iswap++;
    }
  }
  int proc;
  for (int idim = 0; idim < dimension; idim++) {
    noverlap = 0;
    iswap = 2 * idim;
    n = nsendproc[iswap];
    for (i = 0; i < n; i++) {
      proc = sendproc[iswap][i];
      if (proc == me) continue;
      if ((this->*box_touch)(proc, idim, 0)) {
        if (noverlap == maxoverlap) {
          maxoverlap += DELTA_PROCS;
          memory->grow(overlap, maxoverlap, "comm:overlap");
        }
        overlap[noverlap++] = proc;
      }
    }
    noverlap1 = noverlap;
    iswap = 2 * idim + 1;
    n = nsendproc[iswap];
    MPI_Barrier(world);
    for (i = 0; i < n; i++) {
      proc = sendproc[iswap][i];
      if (proc == me) continue;
      if ((this->*box_touch)(proc, idim, 1)) {
        for (j = 0; j < noverlap1; j++)
          if (overlap[j] == proc) break;
        if (j < noverlap1) continue;
        if (noverlap == maxoverlap) {
          maxoverlap += DELTA_PROCS;
          memory->grow(overlap, maxoverlap, "comm:overlap");
        }
        overlap[noverlap++] = proc;
      }
    }
    MPI_Barrier(world);
    if (noverlap > nexchprocmax[idim]) {
      while (nexchprocmax[idim] < noverlap) nexchprocmax[idim] += DELTA_PROCS;
      delete[] exchproc[idim];
      exchproc[idim] = new int[nexchprocmax[idim]];
      delete[] exchnum[idim];
      exchnum[idim] = new int[nexchprocmax[idim]];
    }
    nexchproc[idim] = noverlap;
    for (i = 0; i < noverlap; i++) exchproc[idim][i] = overlap[i];
  }
  if (cutzero) {
    for (i = 0; i < nswap; i++) {
      nsendproc[i] = nrecvproc[i] = sendother[i] = recvother[i] = sendself[i] = 0;
    }
  }
  int nmax = 0;
  for (i = 0; i < nswap; i++) nmax = MAX(nmax, nprocmax[i]);
  for (i = 0; i < dimension; i++) nmax = MAX(nmax, nexchprocmax[i]);
  if (nmax > maxrequest) {
    maxrequest = nmax;
    delete[] requests;
    requests = new MPI_Request[maxrequest];
  }
}
void CommTiled::forward_comm(int )
{
  int i, irecv, n, nsend, nrecv;
  AtomVec *avec = atom->avec;
  double **x = atom->x;
  for (int iswap = 0; iswap < nswap; iswap++) {
    nsend = nsendproc[iswap] - sendself[iswap];
    nrecv = nrecvproc[iswap] - sendself[iswap];
    if (comm_x_only) {
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++)
          MPI_Irecv(x[firstrecv[iswap][i]], size_forward_recv[iswap][i], MPI_DOUBLE,
                    recvproc[iswap][i], 0, world, &requests[i]);
      }
      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++) {
          n = avec->pack_comm(sendnum[iswap][i], sendlist[iswap][i], buf_send, pbc_flag[iswap][i],
                              pbc[iswap][i]);
          MPI_Send(buf_send, n, MPI_DOUBLE, sendproc[iswap][i], 0, world);
        }
      }
      if (sendself[iswap]) {
        avec->pack_comm(sendnum[iswap][nsend], sendlist[iswap][nsend], x[firstrecv[iswap][nrecv]],
                        pbc_flag[iswap][nsend], pbc[iswap][nsend]);
      }
      if (recvother[iswap]) MPI_Waitall(nrecv, requests, MPI_STATUS_IGNORE);
    } else if (ghost_velocity) {
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++)
          MPI_Irecv(&buf_recv[size_forward * forward_recv_offset[iswap][i]],
                    size_forward_recv[iswap][i], MPI_DOUBLE, recvproc[iswap][i], 0, world,
                    &requests[i]);
      }
      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++) {
          n = avec->pack_comm_vel(sendnum[iswap][i], sendlist[iswap][i], buf_send,
                                  pbc_flag[iswap][i], pbc[iswap][i]);
          MPI_Send(buf_send, n, MPI_DOUBLE, sendproc[iswap][i], 0, world);
        }
      }
      if (sendself[iswap]) {
        avec->pack_comm_vel(sendnum[iswap][nsend], sendlist[iswap][nsend], buf_send,
                            pbc_flag[iswap][nsend], pbc[iswap][nsend]);
        avec->unpack_comm_vel(recvnum[iswap][nrecv], firstrecv[iswap][nrecv], buf_send);
      }
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++) {
          MPI_Waitany(nrecv, requests, &irecv, MPI_STATUS_IGNORE);
          avec->unpack_comm_vel(recvnum[iswap][irecv], firstrecv[iswap][irecv],
                                &buf_recv[size_forward * forward_recv_offset[iswap][irecv]]);
        }
      }
    } else {
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++)
          MPI_Irecv(&buf_recv[size_forward * forward_recv_offset[iswap][i]],
                    size_forward_recv[iswap][i], MPI_DOUBLE, recvproc[iswap][i], 0, world,
                    &requests[i]);
      }
      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++) {
          n = avec->pack_comm(sendnum[iswap][i], sendlist[iswap][i], buf_send, pbc_flag[iswap][i],
                              pbc[iswap][i]);
          MPI_Send(buf_send, n, MPI_DOUBLE, sendproc[iswap][i], 0, world);
        }
      }
      if (sendself[iswap]) {
        avec->pack_comm(sendnum[iswap][nsend], sendlist[iswap][nsend], buf_send,
                        pbc_flag[iswap][nsend], pbc[iswap][nsend]);
        avec->unpack_comm(recvnum[iswap][nrecv], firstrecv[iswap][nrecv], buf_send);
      }
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++) {
          MPI_Waitany(nrecv, requests, &irecv, MPI_STATUS_IGNORE);
          avec->unpack_comm(recvnum[iswap][irecv], firstrecv[iswap][irecv],
                            &buf_recv[size_forward * forward_recv_offset[iswap][irecv]]);
        }
      }
    }
  }
}
void CommTiled::reverse_comm()
{
  int i, irecv, n, nsend, nrecv;
  AtomVec *avec = atom->avec;
  double **f = atom->f;
  for (int iswap = nswap - 1; iswap >= 0; iswap--) {
    nsend = nsendproc[iswap] - sendself[iswap];
    nrecv = nrecvproc[iswap] - sendself[iswap];
    if (comm_f_only) {
      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++) {
          MPI_Irecv(&buf_recv[size_reverse * reverse_recv_offset[iswap][i]],
                    size_reverse_recv[iswap][i], MPI_DOUBLE, sendproc[iswap][i], 0, world,
                    &requests[i]);
        }
      }
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++)
          MPI_Send(f[firstrecv[iswap][i]], size_reverse_send[iswap][i], MPI_DOUBLE,
                   recvproc[iswap][i], 0, world);
      }
      if (sendself[iswap]) {
        avec->unpack_reverse(sendnum[iswap][nsend], sendlist[iswap][nsend],
                             f[firstrecv[iswap][nrecv]]);
      }
      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++) {
          MPI_Waitany(nsend, requests, &irecv, MPI_STATUS_IGNORE);
          avec->unpack_reverse(sendnum[iswap][irecv], sendlist[iswap][irecv],
                               &buf_recv[size_reverse * reverse_recv_offset[iswap][irecv]]);
        }
      }
    } else {
      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++)
          MPI_Irecv(&buf_recv[size_reverse * reverse_recv_offset[iswap][i]],
                    size_reverse_recv[iswap][i], MPI_DOUBLE, sendproc[iswap][i], 0, world,
                    &requests[i]);
      }
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++) {
          n = avec->pack_reverse(recvnum[iswap][i], firstrecv[iswap][i], buf_send);
          MPI_Send(buf_send, n, MPI_DOUBLE, recvproc[iswap][i], 0, world);
        }
      }
      if (sendself[iswap]) {
        avec->pack_reverse(recvnum[iswap][nrecv], firstrecv[iswap][nrecv], buf_send);
        avec->unpack_reverse(sendnum[iswap][nsend], sendlist[iswap][nsend], buf_send);
      }
      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++) {
          MPI_Waitany(nsend, requests, &irecv, MPI_STATUS_IGNORE);
          avec->unpack_reverse(sendnum[iswap][irecv], sendlist[iswap][irecv],
                               &buf_recv[size_reverse * reverse_recv_offset[iswap][irecv]]);
        }
      }
    }
  }
}
void CommTiled::exchange()
{
  int i, m, nexch, nsend, nrecv, nlocal, proc, offset;
  double lo, hi, value;
  double **x;
  AtomVec *avec = atom->avec;
  if (map_style != Atom::MAP_NONE) atom->map_clear();
  atom->nghost = 0;
  atom->avec->clear_bonus();
  if (maxexchange_fix_dynamic) {
    int bufextra_old = bufextra;
    init_exchange();
    if (bufextra > bufextra_old) grow_send(maxsend + bufextra, 2);
  }
  if (triclinic == 0) {
    prd = domain->prd;
    boxlo = domain->boxlo;
    boxhi = domain->boxhi;
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    prd = domain->prd_lamda;
    boxlo = domain->boxlo_lamda;
    boxhi = domain->boxhi_lamda;
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }
  dimension = domain->dimension;
  for (int dim = 0; dim < dimension; dim++) {
    x = atom->x;
    lo = sublo[dim];
    hi = subhi[dim];
    nlocal = atom->nlocal;
    i = nsend = 0;
    while (i < nlocal) {
      if (x[i][dim] < lo || x[i][dim] >= hi) {
        if (nsend > maxsend) grow_send(nsend, 1);
        proc = (this->*point_drop)(dim, x[i]);
        if (proc != me) {
          buf_send[nsend++] = proc;
          nsend += avec->pack_exchange(i, &buf_send[nsend]);
        } else {
        }
        avec->copy(nlocal - 1, i, 1);
        nlocal--;
      } else
        i++;
    }
    atom->nlocal = nlocal;
    nexch = nexchproc[dim];
    if (!nexch) continue;
    for (m = 0; m < nexch; m++)
      MPI_Irecv(&exchnum[dim][m], 1, MPI_INT, exchproc[dim][m], 0, world, &requests[m]);
    for (m = 0; m < nexch; m++) MPI_Send(&nsend, 1, MPI_INT, exchproc[dim][m], 0, world);
    MPI_Waitall(nexch, requests, MPI_STATUS_IGNORE);
    nrecv = 0;
    for (m = 0; m < nexch; m++) nrecv += exchnum[dim][m];
    if (nrecv > maxrecv) grow_recv(nrecv);
    offset = 0;
    for (m = 0; m < nexch; m++) {
      MPI_Irecv(&buf_recv[offset], exchnum[dim][m], MPI_DOUBLE, exchproc[dim][m], 0, world,
                &requests[m]);
      offset += exchnum[dim][m];
    }
    for (m = 0; m < nexch; m++) MPI_Send(buf_send, nsend, MPI_DOUBLE, exchproc[dim][m], 0, world);
    MPI_Waitall(nexch, requests, MPI_STATUS_IGNORE);
    m = 0;
    while (m < nrecv) {
      proc = static_cast<int>(buf_recv[m++]);
      if (proc == me) {
        value = buf_recv[m + dim + 1];
        if (value >= lo && value < hi) {
          m += avec->unpack_exchange(&buf_recv[m]);
          continue;
        } else {
        }
      }
      m += static_cast<int>(buf_recv[m]);
    }
  }
  if (atom->firstgroupname) atom->first_reorder();
}
void CommTiled::borders()
{
  int i, m, n, nlast, nsend, nrecv, ngroup, nprior, ncount, ncountall;
  double xlo, xhi, ylo, yhi, zlo, zhi;
  double *bbox;
  double **x;
  AtomVec *avec = atom->avec;
  if (mode == Comm::MULTI) neighbor->build_collection(0);
  smaxone = smaxall = 0;
  rmaxone = rmaxall = 0;
  for (int iswap = 0; iswap < nswap; iswap++) {
    x = atom->x;
    if (iswap % 2 == 0) nlast = atom->nlocal + atom->nghost;
    ncountall = 0;
    for (m = 0; m < nsendproc[iswap]; m++) {
      if (mode == Comm::SINGLE) {
        bbox = sendbox[iswap][m];
        xlo = bbox[0];
        ylo = bbox[1];
        zlo = bbox[2];
        xhi = bbox[3];
        yhi = bbox[4];
        zhi = bbox[5];
        ncount = 0;
        if (!bordergroup) {
          for (i = 0; i < nlast; i++) {
            if (x[i][0] >= xlo && x[i][0] < xhi && x[i][1] >= ylo && x[i][1] < yhi &&
                x[i][2] >= zlo && x[i][2] < zhi) {
              if (ncount == maxsendlist[iswap][m]) grow_list(iswap, m, ncount);
              sendlist[iswap][m][ncount++] = i;
            }
          }
        } else {
          ngroup = atom->nfirst;
          for (i = 0; i < ngroup; i++) {
            if (x[i][0] >= xlo && x[i][0] < xhi && x[i][1] >= ylo && x[i][1] < yhi &&
                x[i][2] >= zlo && x[i][2] < zhi) {
              if (ncount == maxsendlist[iswap][m]) grow_list(iswap, m, ncount);
              sendlist[iswap][m][ncount++] = i;
            }
          }
          for (i = atom->nlocal; i < nlast; i++) {
            if (x[i][0] >= xlo && x[i][0] < xhi && x[i][1] >= ylo && x[i][1] < yhi &&
                x[i][2] >= zlo && x[i][2] < zhi) {
              if (ncount == maxsendlist[iswap][m]) grow_list(iswap, m, ncount);
              sendlist[iswap][m][ncount++] = i;
            }
          }
        }
        sendnum[iswap][m] = ncount;
        smaxone = MAX(smaxone, ncount);
        ncountall += ncount;
      } else if (mode == Comm::MULTI) {
        int *collection = neighbor->collection;
        int icollection;
        ncount = 0;
        if (!bordergroup) {
          for (i = 0; i < nlast; i++) {
            icollection = collection[i];
            bbox = sendbox_multi[iswap][m][icollection];
            xlo = bbox[0];
            ylo = bbox[1];
            zlo = bbox[2];
            xhi = bbox[3];
            yhi = bbox[4];
            zhi = bbox[5];
            if (x[i][0] >= xlo && x[i][0] < xhi && x[i][1] >= ylo && x[i][1] < yhi &&
                x[i][2] >= zlo && x[i][2] < zhi) {
              if (ncount == maxsendlist[iswap][m]) grow_list(iswap, m, ncount);
              sendlist[iswap][m][ncount++] = i;
            }
          }
        } else {
          ngroup = atom->nfirst;
          for (i = 0; i < ngroup; i++) {
            icollection = collection[i];
            bbox = sendbox_multi[iswap][m][icollection];
            xlo = bbox[0];
            ylo = bbox[1];
            zlo = bbox[2];
            xhi = bbox[3];
            yhi = bbox[4];
            zhi = bbox[5];
            if (x[i][0] >= xlo && x[i][0] < xhi && x[i][1] >= ylo && x[i][1] < yhi &&
                x[i][2] >= zlo && x[i][2] < zhi) {
              if (ncount == maxsendlist[iswap][m]) grow_list(iswap, m, ncount);
              sendlist[iswap][m][ncount++] = i;
            }
          }
          for (i = atom->nlocal; i < nlast; i++) {
            icollection = collection[i];
            bbox = sendbox_multi[iswap][m][icollection];
            xlo = bbox[0];
            ylo = bbox[1];
            zlo = bbox[2];
            xhi = bbox[3];
            yhi = bbox[4];
            zhi = bbox[5];
            if (x[i][0] >= xlo && x[i][0] < xhi && x[i][1] >= ylo && x[i][1] < yhi &&
                x[i][2] >= zlo && x[i][2] < zhi) {
              if (ncount == maxsendlist[iswap][m]) grow_list(iswap, m, ncount);
              sendlist[iswap][m][ncount++] = i;
            }
          }
        }
        sendnum[iswap][m] = ncount;
        smaxone = MAX(smaxone, ncount);
        ncountall += ncount;
      } else {
        int *type = atom->type;
        int itype;
        ncount = 0;
        if (!bordergroup) {
          for (i = 0; i < nlast; i++) {
            itype = type[i];
            bbox = sendbox_multiold[iswap][m][itype];
            xlo = bbox[0];
            ylo = bbox[1];
            zlo = bbox[2];
            xhi = bbox[3];
            yhi = bbox[4];
            zhi = bbox[5];
            if (x[i][0] >= xlo && x[i][0] < xhi && x[i][1] >= ylo && x[i][1] < yhi &&
                x[i][2] >= zlo && x[i][2] < zhi) {
              if (ncount == maxsendlist[iswap][m]) grow_list(iswap, m, ncount);
              sendlist[iswap][m][ncount++] = i;
            }
          }
        } else {
          ngroup = atom->nfirst;
          for (i = 0; i < ngroup; i++) {
            itype = type[i];
            bbox = sendbox_multiold[iswap][m][itype];
            xlo = bbox[0];
            ylo = bbox[1];
            zlo = bbox[2];
            xhi = bbox[3];
            yhi = bbox[4];
            zhi = bbox[5];
            if (x[i][0] >= xlo && x[i][0] < xhi && x[i][1] >= ylo && x[i][1] < yhi &&
                x[i][2] >= zlo && x[i][2] < zhi) {
              if (ncount == maxsendlist[iswap][m]) grow_list(iswap, m, ncount);
              sendlist[iswap][m][ncount++] = i;
            }
          }
          for (i = atom->nlocal; i < nlast; i++) {
            itype = type[i];
            bbox = sendbox_multiold[iswap][m][itype];
            xlo = bbox[0];
            ylo = bbox[1];
            zlo = bbox[2];
            xhi = bbox[3];
            yhi = bbox[4];
            zhi = bbox[5];
            if (x[i][0] >= xlo && x[i][0] < xhi && x[i][1] >= ylo && x[i][1] < yhi &&
                x[i][2] >= zlo && x[i][2] < zhi) {
              if (ncount == maxsendlist[iswap][m]) grow_list(iswap, m, ncount);
              sendlist[iswap][m][ncount++] = i;
            }
          }
        }
        sendnum[iswap][m] = ncount;
        smaxone = MAX(smaxone, ncount);
        ncountall += ncount;
      }
    }
    smaxall = MAX(smaxall, ncountall);
    nsend = nsendproc[iswap] - sendself[iswap];
    nrecv = nrecvproc[iswap] - sendself[iswap];
    if (recvother[iswap])
      for (m = 0; m < nrecv; m++)
        MPI_Irecv(&recvnum[iswap][m], 1, MPI_INT, recvproc[iswap][m], 0, world, &requests[m]);
    if (sendother[iswap])
      for (m = 0; m < nsend; m++)
        MPI_Send(&sendnum[iswap][m], 1, MPI_INT, sendproc[iswap][m], 0, world);
    if (sendself[iswap]) recvnum[iswap][nrecv] = sendnum[iswap][nsend];
    if (recvother[iswap]) MPI_Waitall(nrecv, requests, MPI_STATUS_IGNORE);
    for (m = 0; m < nsendproc[iswap]; m++) {
      size_reverse_recv[iswap][m] = sendnum[iswap][m] * size_reverse;
      if (m == 0)
        reverse_recv_offset[iswap][0] = 0;
      else
        reverse_recv_offset[iswap][m] = reverse_recv_offset[iswap][m - 1] + sendnum[iswap][m - 1];
    }
    ncountall = 0;
    for (m = 0; m < nrecvproc[iswap]; m++) {
      ncount = recvnum[iswap][m];
      rmaxone = MAX(rmaxone, ncount);
      ncountall += ncount;
      size_forward_recv[iswap][m] = ncount * size_forward;
      size_reverse_send[iswap][m] = ncount * size_reverse;
      if (m == 0) {
        firstrecv[iswap][0] = atom->nlocal + atom->nghost;
        forward_recv_offset[iswap][0] = 0;
      } else {
        firstrecv[iswap][m] = firstrecv[iswap][m - 1] + recvnum[iswap][m - 1];
        forward_recv_offset[iswap][m] = forward_recv_offset[iswap][m - 1] + recvnum[iswap][m - 1];
      }
    }
    rmaxall = MAX(rmaxall, ncountall);
    if (smaxone * size_border > maxsend) grow_send(smaxone * size_border, 0);
    if (rmaxall * size_border > maxrecv) grow_recv(rmaxall * size_border);
    if (ghost_velocity) {
      if (recvother[iswap]) {
        for (m = 0; m < nrecv; m++)
          MPI_Irecv(&buf_recv[size_border * forward_recv_offset[iswap][m]],
                    recvnum[iswap][m] * size_border, MPI_DOUBLE, recvproc[iswap][m], 0, world,
                    &requests[m]);
      }
      if (sendother[iswap]) {
        for (m = 0; m < nsend; m++) {
          n = avec->pack_border_vel(sendnum[iswap][m], sendlist[iswap][m], buf_send,
                                    pbc_flag[iswap][m], pbc[iswap][m]);
          MPI_Send(buf_send, n, MPI_DOUBLE, sendproc[iswap][m], 0, world);
        }
      }
      if (sendself[iswap]) {
        avec->pack_border_vel(sendnum[iswap][nsend], sendlist[iswap][nsend], buf_send,
                              pbc_flag[iswap][nsend], pbc[iswap][nsend]);
        avec->unpack_border_vel(recvnum[iswap][nrecv], firstrecv[iswap][nrecv], buf_send);
      }
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++) {
          MPI_Waitany(nrecv, requests, &m, MPI_STATUS_IGNORE);
          avec->unpack_border_vel(recvnum[iswap][m], firstrecv[iswap][m],
                                  &buf_recv[size_border * forward_recv_offset[iswap][m]]);
        }
      }
    } else {
      if (recvother[iswap]) {
        for (m = 0; m < nrecv; m++)
          MPI_Irecv(&buf_recv[size_border * forward_recv_offset[iswap][m]],
                    recvnum[iswap][m] * size_border, MPI_DOUBLE, recvproc[iswap][m], 0, world,
                    &requests[m]);
      }
      if (sendother[iswap]) {
        for (m = 0; m < nsend; m++) {
          n = avec->pack_border(sendnum[iswap][m], sendlist[iswap][m], buf_send, pbc_flag[iswap][m],
                                pbc[iswap][m]);
          MPI_Send(buf_send, n, MPI_DOUBLE, sendproc[iswap][m], 0, world);
        }
      }
      if (sendself[iswap]) {
        avec->pack_border(sendnum[iswap][nsend], sendlist[iswap][nsend], buf_send,
                          pbc_flag[iswap][nsend], pbc[iswap][nsend]);
        avec->unpack_border(recvnum[iswap][nrecv], firstrecv[iswap][nrecv], buf_send);
      }
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++) {
          MPI_Waitany(nrecv, requests, &m, MPI_STATUS_IGNORE);
          avec->unpack_border(recvnum[iswap][m], firstrecv[iswap][m],
                              &buf_recv[size_border * forward_recv_offset[iswap][m]]);
        }
      }
    }
    n = nrecvproc[iswap];
    if (n) {
      nprior = atom->nghost + atom->nlocal;
      atom->nghost += forward_recv_offset[iswap][n - 1] + recvnum[iswap][n - 1];
      if (neighbor->style == Neighbor::MULTI) neighbor->build_collection(nprior);
    }
  }
  int max = MAX(maxforward * smaxone, maxreverse * rmaxone);
  if (max > maxsend) grow_send(max, 0);
  max = MAX(maxforward * rmaxall, maxreverse * smaxall);
  if (max > maxrecv) grow_recv(max);
  if (map_style != Atom::MAP_NONE) atom->map_set();
}
void CommTiled::forward_comm(Pair *pair)
{
  int i, irecv, n, nsend, nrecv;
  int nsize = pair->comm_forward;
  for (int iswap = 0; iswap < nswap; iswap++) {
    nsend = nsendproc[iswap] - sendself[iswap];
    nrecv = nrecvproc[iswap] - sendself[iswap];
    if (recvother[iswap]) {
      for (i = 0; i < nrecv; i++)
        MPI_Irecv(&buf_recv[nsize * forward_recv_offset[iswap][i]], nsize * recvnum[iswap][i],
                  MPI_DOUBLE, recvproc[iswap][i], 0, world, &requests[i]);
    }
    if (sendother[iswap]) {
      for (i = 0; i < nsend; i++) {
        n = pair->pack_forward_comm(sendnum[iswap][i], sendlist[iswap][i], buf_send,
                                    pbc_flag[iswap][i], pbc[iswap][i]);
        MPI_Send(buf_send, n, MPI_DOUBLE, sendproc[iswap][i], 0, world);
      }
    }
    if (sendself[iswap]) {
      pair->pack_forward_comm(sendnum[iswap][nsend], sendlist[iswap][nsend], buf_send,
                              pbc_flag[iswap][nsend], pbc[iswap][nsend]);
      pair->unpack_forward_comm(recvnum[iswap][nrecv], firstrecv[iswap][nrecv], buf_send);
    }
    if (recvother[iswap]) {
      for (i = 0; i < nrecv; i++) {
        MPI_Waitany(nrecv, requests, &irecv, MPI_STATUS_IGNORE);
        pair->unpack_forward_comm(recvnum[iswap][irecv], firstrecv[iswap][irecv],
                                  &buf_recv[nsize * forward_recv_offset[iswap][irecv]]);
      }
    }
  }
}
void CommTiled::reverse_comm(Pair *pair)
{
  int i, irecv, n, nsend, nrecv;
  int nsize = MAX(pair->comm_reverse, pair->comm_reverse_off);
  for (int iswap = nswap - 1; iswap >= 0; iswap--) {
    nsend = nsendproc[iswap] - sendself[iswap];
    nrecv = nrecvproc[iswap] - sendself[iswap];
    if (sendother[iswap]) {
      for (i = 0; i < nsend; i++)
        MPI_Irecv(&buf_recv[nsize * reverse_recv_offset[iswap][i]], nsize * sendnum[iswap][i],
                  MPI_DOUBLE, sendproc[iswap][i], 0, world, &requests[i]);
    }
    if (recvother[iswap]) {
      for (i = 0; i < nrecv; i++) {
        n = pair->pack_reverse_comm(recvnum[iswap][i], firstrecv[iswap][i], buf_send);
        MPI_Send(buf_send, n, MPI_DOUBLE, recvproc[iswap][i], 0, world);
      }
    }
    if (sendself[iswap]) {
      pair->pack_reverse_comm(recvnum[iswap][nrecv], firstrecv[iswap][nrecv], buf_send);
      pair->unpack_reverse_comm(sendnum[iswap][nsend], sendlist[iswap][nsend], buf_send);
    }
    if (sendother[iswap]) {
      for (i = 0; i < nsend; i++) {
        MPI_Waitany(nsend, requests, &irecv, MPI_STATUS_IGNORE);
        pair->unpack_reverse_comm(sendnum[iswap][irecv], sendlist[iswap][irecv],
                                  &buf_recv[nsize * reverse_recv_offset[iswap][irecv]]);
      }
    }
  }
}
void CommTiled::forward_comm(Fix *fix, int size)
{
  int i, irecv, n, nsize, nsend, nrecv;
  if (size)
    nsize = size;
  else
    nsize = fix->comm_forward;
  for (int iswap = 0; iswap < nswap; iswap++) {
    nsend = nsendproc[iswap] - sendself[iswap];
    nrecv = nrecvproc[iswap] - sendself[iswap];
    if (recvother[iswap]) {
      for (i = 0; i < nrecv; i++)
        MPI_Irecv(&buf_recv[nsize * forward_recv_offset[iswap][i]], nsize * recvnum[iswap][i],
                  MPI_DOUBLE, recvproc[iswap][i], 0, world, &requests[i]);
    }
    if (sendother[iswap]) {
      for (i = 0; i < nsend; i++) {
        n = fix->pack_forward_comm(sendnum[iswap][i], sendlist[iswap][i], buf_send,
                                   pbc_flag[iswap][i], pbc[iswap][i]);
        MPI_Send(buf_send, n, MPI_DOUBLE, sendproc[iswap][i], 0, world);
      }
    }
    if (sendself[iswap]) {
      fix->pack_forward_comm(sendnum[iswap][nsend], sendlist[iswap][nsend], buf_send,
                             pbc_flag[iswap][nsend], pbc[iswap][nsend]);
      fix->unpack_forward_comm(recvnum[iswap][nrecv], firstrecv[iswap][nrecv], buf_send);
    }
    if (recvother[iswap]) {
      for (i = 0; i < nrecv; i++) {
        MPI_Waitany(nrecv, requests, &irecv, MPI_STATUS_IGNORE);
        fix->unpack_forward_comm(recvnum[iswap][irecv], firstrecv[iswap][irecv],
                                 &buf_recv[nsize * forward_recv_offset[iswap][irecv]]);
      }
    }
  }
}
void CommTiled::reverse_comm(Fix *fix, int size)
{
  int i, irecv, n, nsize, nsend, nrecv;
  if (size)
    nsize = size;
  else
    nsize = fix->comm_reverse;
  for (int iswap = nswap - 1; iswap >= 0; iswap--) {
    nsend = nsendproc[iswap] - sendself[iswap];
    nrecv = nrecvproc[iswap] - sendself[iswap];
    if (sendother[iswap]) {
      for (i = 0; i < nsend; i++)
        MPI_Irecv(&buf_recv[nsize * reverse_recv_offset[iswap][i]], nsize * sendnum[iswap][i],
                  MPI_DOUBLE, sendproc[iswap][i], 0, world, &requests[i]);
    }
    if (recvother[iswap]) {
      for (i = 0; i < nrecv; i++) {
        n = fix->pack_reverse_comm(recvnum[iswap][i], firstrecv[iswap][i], buf_send);
        MPI_Send(buf_send, n, MPI_DOUBLE, recvproc[iswap][i], 0, world);
      }
    }
    if (sendself[iswap]) {
      fix->pack_reverse_comm(recvnum[iswap][nrecv], firstrecv[iswap][nrecv], buf_send);
      fix->unpack_reverse_comm(sendnum[iswap][nsend], sendlist[iswap][nsend], buf_send);
    }
    if (sendother[iswap]) {
      for (i = 0; i < nsend; i++) {
        MPI_Waitany(nsend, requests, &irecv, MPI_STATUS_IGNORE);
        fix->unpack_reverse_comm(sendnum[iswap][irecv], sendlist[iswap][irecv],
                                 &buf_recv[nsize * reverse_recv_offset[iswap][irecv]]);
      }
    }
  }
}
void CommTiled::reverse_comm_variable(Fix * )
{
  error->all(FLERR, "Reverse comm fix variable not yet supported by CommTiled");
}
void CommTiled::forward_comm(Compute *compute)
{
  int i, irecv, n, nsend, nrecv;
  int nsize = compute->comm_forward;
  for (int iswap = 0; iswap < nswap; iswap++) {
    nsend = nsendproc[iswap] - sendself[iswap];
    nrecv = nrecvproc[iswap] - sendself[iswap];
    if (recvother[iswap]) {
      for (i = 0; i < nrecv; i++)
        MPI_Irecv(&buf_recv[nsize * forward_recv_offset[iswap][i]], nsize * recvnum[iswap][i],
                  MPI_DOUBLE, recvproc[iswap][i], 0, world, &requests[i]);
    }
    if (sendother[iswap]) {
      for (i = 0; i < nsend; i++) {
        n = compute->pack_forward_comm(sendnum[iswap][i], sendlist[iswap][i], buf_send,
                                       pbc_flag[iswap][i], pbc[iswap][i]);
        MPI_Send(buf_send, n, MPI_DOUBLE, sendproc[iswap][i], 0, world);
      }
    }
    if (sendself[iswap]) {
      compute->pack_forward_comm(sendnum[iswap][nsend], sendlist[iswap][nsend], buf_send,
                                 pbc_flag[iswap][nsend], pbc[iswap][nsend]);
      compute->unpack_forward_comm(recvnum[iswap][nrecv], firstrecv[iswap][nrecv], buf_send);
    }
    if (recvother[iswap]) {
      for (i = 0; i < nrecv; i++) {
        MPI_Waitany(nrecv, requests, &irecv, MPI_STATUS_IGNORE);
        compute->unpack_forward_comm(recvnum[iswap][irecv], firstrecv[iswap][irecv],
                                     &buf_recv[nsize * forward_recv_offset[iswap][irecv]]);
      }
    }
  }
}
void CommTiled::reverse_comm(Compute *compute)
{
  int i, irecv, n, nsend, nrecv;
  int nsize = compute->comm_reverse;
  for (int iswap = nswap - 1; iswap >= 0; iswap--) {
    nsend = nsendproc[iswap] - sendself[iswap];
    nrecv = nrecvproc[iswap] - sendself[iswap];
    if (sendother[iswap]) {
      for (i = 0; i < nsend; i++)
        MPI_Irecv(&buf_recv[nsize * reverse_recv_offset[iswap][i]], nsize * sendnum[iswap][i],
                  MPI_DOUBLE, sendproc[iswap][i], 0, world, &requests[i]);
    }
    if (recvother[iswap]) {
      for (i = 0; i < nrecv; i++) {
        n = compute->pack_reverse_comm(recvnum[iswap][i], firstrecv[iswap][i], buf_send);
        MPI_Send(buf_send, n, MPI_DOUBLE, recvproc[iswap][i], 0, world);
      }
    }
    if (sendself[iswap]) {
      compute->pack_reverse_comm(recvnum[iswap][nrecv], firstrecv[iswap][nrecv], buf_send);
      compute->unpack_reverse_comm(sendnum[iswap][nsend], sendlist[iswap][nsend], buf_send);
    }
    if (sendother[iswap]) {
      for (i = 0; i < nsend; i++) {
        MPI_Waitany(nsend, requests, &irecv, MPI_STATUS_IGNORE);
        compute->unpack_reverse_comm(sendnum[iswap][irecv], sendlist[iswap][irecv],
                                     &buf_recv[nsize * reverse_recv_offset[iswap][irecv]]);
      }
    }
  }
}
void CommTiled::forward_comm_array(int nsize, double **array)
{
  int i, j, k, m, iatom, last, irecv, nsend, nrecv;
  if (nsize > maxforward) {
    maxforward = nsize;
    if (maxforward * smaxone > maxsend) grow_send(maxforward * smaxone, 0);
    if (maxforward * rmaxall > maxrecv) grow_recv(maxforward * rmaxall);
  }
  for (int iswap = 0; iswap < nswap; iswap++) {
    nsend = nsendproc[iswap] - sendself[iswap];
    nrecv = nrecvproc[iswap] - sendself[iswap];
    MPI_Barrier(world);
    if (recvother[iswap]) {
      for (i = 0; i < nrecv; i++)
        MPI_Irecv(&buf_recv[nsize * forward_recv_offset[iswap][i]], nsize * recvnum[iswap][i],
                  MPI_DOUBLE, recvproc[iswap][i], 0, world, &requests[i]);
    }
    if (sendother[iswap]) {
      for (i = 0; i < nsend; i++) {
        m = 0;
        for (iatom = 0; iatom < sendnum[iswap][i]; iatom++) {
          j = sendlist[iswap][i][iatom];
          for (k = 0; k < nsize; k++) buf_send[m++] = array[j][k];
        }
        MPI_Send(buf_send, nsize * sendnum[iswap][i], MPI_DOUBLE, sendproc[iswap][i], 0, world);
      }
    }
    if (sendself[iswap]) {
      m = 0;
      for (iatom = 0; iatom < sendnum[iswap][nsend]; iatom++) {
        j = sendlist[iswap][nsend][iatom];
        for (k = 0; k < nsize; k++) buf_send[m++] = array[j][k];
      }
      m = 0;
      last = firstrecv[iswap][nrecv] + recvnum[iswap][nrecv];
      for (iatom = firstrecv[iswap][nrecv]; iatom < last; iatom++)
        for (k = 0; k < nsize; k++) array[iatom][k] = buf_send[m++];
    }
    if (recvother[iswap]) {
      for (i = 0; i < nrecv; i++) {
        MPI_Waitany(nrecv, requests, &irecv, MPI_STATUS_IGNORE);
        m = nsize * forward_recv_offset[iswap][irecv];
        last = firstrecv[iswap][irecv] + recvnum[iswap][irecv];
        for (iatom = firstrecv[iswap][irecv]; iatom < last; iatom++)
          for (k = 0; k < nsize; k++) array[iatom][k] = buf_recv[m++];
      }
    }
  }
}
struct Box {
  int rank;
  double lo[3];
  double hi[3];
};
static void dump(int me, int nprocs, int nneigh, int *neighbors, struct Box *boxes)
{
  char path[FILENAME_MAX];
  int i, j, k, Neighbor, ineigh;
  const int shift[8][3] = {
      {0, 0, 0}, {0, 0, 1}, {0, 1, 1}, {0, 1, 0}, {1, 0, 0}, {1, 0, 1}, {1, 1, 1}, {1, 1, 0},
  };
  const int faces[6][4] = {
      {0, 3, 2, 1}, {0, 1, 5, 4}, {0, 4, 7, 3}, {1, 5, 6, 2}, {3, 7, 6, 2}, {4, 5, 6, 7},
  };
  FILE *file;
  snprintf(path, sizeof path, "box.%05d.off", me);
  if ((file = fopen(path, "w")) == NULL) {
    fprintf(stderr, "fopen failed for '%s'\n", path);
    exit(1);
  }
  fprintf(file, "OFF\n");
  fprintf(file, "%d %d 0\n", 8 * nprocs, 6 * nprocs);
  for (k = 0; k < nprocs; k++) {
    for (j = 0; j < 8; j++) {
      for (i = 0; i < 3; i++)
        fprintf(file, "%.16e ", shift[j][i] ? boxes[k].hi[i] : boxes[k].lo[i]);
      fprintf(file, "\n");
    }
  }
  for (k = 0; k < nprocs; k++) {
    for (j = 0; j < 6; j++) {
      fprintf(file, "4");
      for (i = 0; i < 4; i++) fprintf(file, " %d", 8 * k + faces[j][i]);
      if (boxes[k].rank == me)
        fprintf(file, " 1 0 0");
      else {
        Neighbor = 0;
        for (ineigh = 0; ineigh < nneigh; ineigh++)
          if (boxes[k].rank == neighbors[ineigh]) {
            fprintf(file, " 0 1 0");
            Neighbor = 1;
            break;
          }
        if (!Neighbor) fprintf(file, " 0.75 0.75 0.75");
      }
      fprintf(file, "\n");
    }
  }
  if (fclose(file) != 0) {
    fprintf(stderr, "fclose failed for %s\n", path);
    exit(1);
  }
}
int CommTiled::exchange_variable(int n, double *inbuf, double *&outbuf)
{
  Box box, *boxes;
  int d, k, s, mneigh;
  double lo, hi, u, v, dist, dist_sq, skin = 12;
  int ineigh;
  static int Ready = 0;
  static int *neighbors;
  static int nneigh;
  if (Ready == 0) {
    if ((boxes = (struct Box *) malloc(nprocs * sizeof *boxes)) == NULL)
      error->one(FLERR, "malloc faild\n");
    box.lo[0] = sublo[0];
    box.lo[1] = sublo[1];
    box.lo[2] = sublo[2];
    box.hi[0] = subhi[0];
    box.hi[1] = subhi[1];
    box.hi[2] = subhi[2];
    box.rank = me;
    MPI_Allgather(&box, sizeof box, MPI_BYTE, boxes, sizeof box, MPI_BYTE, world);
    neighbors = NULL;
    nneigh = 0;
    mneigh = 0;
    if (domain->triclinic) error->one(FLERR, "does not work for triclinic\n");
    for (k = 0; k < nprocs; k++)
      if (boxes[k].rank != me) {
        dist_sq = 0;
        for (d = 0; d < 3; d++) {
          dist = DBL_MAX;
          for (s = -1; s <= 1; s++)
            if (s == 0 || domain->periodicity[d]) {
              lo = sublo[d] + s * domain->prd[d];
              hi = subhi[d] + s * domain->prd[d];
              if (lo < boxes[k].lo[d]) {
                u = hi - lo;
                v = boxes[k].lo[d] - lo;
              } else {
                u = boxes[k].hi[d] - boxes[k].lo[d];
                v = lo - boxes[k].lo[d];
              }
              if (v > u) {
                if (v - u < dist) dist = v - u;
              } else
                dist = 0;
            }
          dist_sq += dist * dist;
        }
        if (dist_sq <= skin * skin) {
          if (nneigh >= mneigh) {
            mneigh = 2 * mneigh + 1;
            if ((neighbors = (int *) realloc(neighbors, mneigh * sizeof *neighbors)) == NULL)
              error->one(FLERR, "realloc faild\n");
          }
          neighbors[nneigh] = boxes[k].rank;
          nneigh++;
        }
      }
    free(boxes);
    Ready = 1;
  }
  int *recived, flag, offset, nmsg;
  MPI_Status status;
  MPI_Request *req;
  if ((req = (MPI_Request *) malloc(nneigh * sizeof *req)) == NULL)
    error->one(FLERR, "malloc faild\n");
  for (ineigh = 0; ineigh < nneigh; ineigh++)
    MPI_Isend(inbuf, n, MPI_DOUBLE, neighbors[ineigh], 0, world, &req[ineigh]);
  if (n > maxrecv) grow_recv(n);
  if (inbuf) memcpy(buf_recv, inbuf, n * sizeof *inbuf);
  recived = NULL;
  if ((recived = (int *) calloc(nneigh, sizeof *recived)) == NULL)
    error->one(FLERR, "calloc failed");
  nmsg = 0;
  for (;;) {
    for (ineigh = 0; ineigh < nneigh; ineigh++)
      if (!recived[ineigh]) {
        MPI_Iprobe(neighbors[ineigh], 0, world, &flag, &status);
        if (flag) {
          recived[ineigh] = 1;
          MPI_Get_count(&status, MPI_DOUBLE, &k);
          offset = n;
          n += k;
          if (n > maxrecv) {
            maxrecv = static_cast<int>(BUFFACTOR * n);
            memory->grow(buf_recv, maxrecv + bufextra, "comm:buf_recv");
          }
          MPI_Recv(&buf_recv[offset], k, MPI_DOUBLE, neighbors[ineigh], 0, world,
                   MPI_STATUS_IGNORE);
          nmsg++;
          if (nmsg == nneigh) goto end;
        }
      }
  }
end:
  MPI_Waitall(nneigh, req, MPI_STATUS_IGNORE);
  free(recived);
  free(req);
  outbuf = buf_recv;
  return n;
}
int CommTiled::exchange_variable_all2all(int nsend, double *inbuf, double *&outbuf)
{
  error->one(FLERR, "CommBrick::exchange_variable_all2all() is not implimented");
}
void CommTiled::box_drop_brick(int idim, double *lo, double *hi, int &indexme)
{
  int dir;
  int index = -1;
  if (hi[idim] == sublo[idim]) {
    index = myloc[idim] - 1;
    dir = -1;
  } else if (lo[idim] == subhi[idim]) {
    index = myloc[idim] + 1;
    dir = 1;
  } else if (hi[idim] == boxhi[idim]) {
    index = procgrid[idim] - 1;
    dir = -1;
  } else if (lo[idim] == boxlo[idim]) {
    index = 0;
    dir = 1;
  } else
    error->one(FLERR, "Comm tiled mis-match in box drop brick");
  int other1, other2, proc;
  double lower, upper;
  double *split;
  if (idim == 0) {
    other1 = myloc[1];
    other2 = myloc[2];
    split = xsplit;
  } else if (idim == 1) {
    other1 = myloc[0];
    other2 = myloc[2];
    split = ysplit;
  } else {
    other1 = myloc[0];
    other2 = myloc[1];
    split = zsplit;
  }
  if (index < 0 || index > procgrid[idim])
    error->one(FLERR, "Comm tiled invalid index in box drop brick");
  while (true) {
    lower = boxlo[idim] + prd[idim] * split[index];
    if (index < procgrid[idim] - 1)
      upper = boxlo[idim] + prd[idim] * split[index + 1];
    else
      upper = boxhi[idim];
    if (lower >= hi[idim] || upper <= lo[idim]) break;
    if (idim == 0)
      proc = grid2proc[index][other1][other2];
    else if (idim == 1)
      proc = grid2proc[other1][index][other2];
    else
      proc = grid2proc[other1][other2][index];
    if (noverlap == maxoverlap) {
      maxoverlap += DELTA_PROCS;
      memory->grow(overlap, maxoverlap, "comm:overlap");
    }
    if (proc == me) indexme = noverlap;
    overlap[noverlap++] = proc;
    index += dir;
    if (index < 0 || index >= procgrid[idim]) break;
  }
}
void CommTiled::box_drop_tiled(int , double *lo, double *hi, int &indexme)
{
  box_drop_tiled_recurse(lo, hi, 0, nprocs - 1, indexme);
}
void CommTiled::box_drop_tiled_recurse(double *lo, double *hi, int proclower, int procupper,
                                       int &indexme)
{
  if (proclower == procupper) {
    if (noverlap == maxoverlap) {
      maxoverlap += DELTA_PROCS;
      memory->grow(overlap, maxoverlap, "comm:overlap");
    }
    if (proclower == me) indexme = noverlap;
    overlap[noverlap++] = proclower;
    return;
  }
  int procmid = proclower + (procupper - proclower) / 2 + 1;
  int idim = rcbinfo[procmid].dim;
  double cut = boxlo[idim] + prd[idim] * rcbinfo[procmid].cutfrac;
  if (lo[idim] < cut) box_drop_tiled_recurse(lo, hi, proclower, procmid - 1, indexme);
  if (hi[idim] > cut) box_drop_tiled_recurse(lo, hi, procmid, procupper, indexme);
}
void CommTiled::box_other_brick(int idim, int idir, int proc, double *lo, double *hi)
{
  lo[0] = sublo[0];
  lo[1] = sublo[1];
  lo[2] = sublo[2];
  hi[0] = subhi[0];
  hi[1] = subhi[1];
  hi[2] = subhi[2];
  int other1, other2, oproc;
  double *split;
  if (idim == 0) {
    other1 = myloc[1];
    other2 = myloc[2];
    split = xsplit;
  } else if (idim == 1) {
    other1 = myloc[0];
    other2 = myloc[2];
    split = ysplit;
  } else {
    other1 = myloc[0];
    other2 = myloc[1];
    split = zsplit;
  }
  int dir = -1;
  if (idir) dir = 1;
  int index = myloc[idim];
  int n = procgrid[idim];
  for (int i = 0; i < n; i++) {
    index += dir;
    if (index < 0)
      index = n - 1;
    else if (index >= n)
      index = 0;
    if (idim == 0)
      oproc = grid2proc[index][other1][other2];
    else if (idim == 1)
      oproc = grid2proc[other1][index][other2];
    else
      oproc = grid2proc[other1][other2][index];
    if (proc == oproc) {
      lo[idim] = boxlo[idim] + prd[idim] * split[index];
      if (split[index + 1] < 1.0)
        hi[idim] = boxlo[idim] + prd[idim] * split[index + 1];
      else
        hi[idim] = boxhi[idim];
      return;
    }
  }
}
void CommTiled::box_other_tiled(int , int , int proc, double *lo, double *hi)
{
  double(*split)[2] = rcbinfo[proc].mysplit;
  lo[0] = boxlo[0] + prd[0] * split[0][0];
  if (split[0][1] < 1.0)
    hi[0] = boxlo[0] + prd[0] * split[0][1];
  else
    hi[0] = boxhi[0];
  lo[1] = boxlo[1] + prd[1] * split[1][0];
  if (split[1][1] < 1.0)
    hi[1] = boxlo[1] + prd[1] * split[1][1];
  else
    hi[1] = boxhi[1];
  lo[2] = boxlo[2] + prd[2] * split[2][0];
  if (split[2][1] < 1.0)
    hi[2] = boxlo[2] + prd[2] * split[2][1];
  else
    hi[2] = boxhi[2];
}
int CommTiled::box_touch_brick(int proc, int idim, int idir)
{
  if (procneigh[idim][idir] == proc) return 1;
  return 0;
}
int CommTiled::box_touch_tiled(int proc, int idim, int idir)
{
  if (idir == 0) {
    if (rcbinfo[proc].mysplit[idim][1] == rcbinfo[me].mysplit[idim][0])
      return 1;
    else if (rcbinfo[proc].mysplit[idim][1] == 1.0 && rcbinfo[me].mysplit[idim][0] == 0.0)
      return 1;
  } else {
    if (rcbinfo[proc].mysplit[idim][0] == rcbinfo[me].mysplit[idim][1])
      return 1;
    else if (rcbinfo[proc].mysplit[idim][0] == 0.0 && rcbinfo[me].mysplit[idim][1] == 1.0)
      return 1;
  }
  return 0;
}
int CommTiled::point_drop_brick(int idim, double *x)
{
  if (closer_subbox_edge(idim, x)) return procneigh[idim][1];
  return procneigh[idim][0];
}
int CommTiled::point_drop_tiled(int idim, double *x)
{
  double xnew[3];
  xnew[0] = x[0];
  xnew[1] = x[1];
  xnew[2] = x[2];
  if (idim == 0) {
    if (xnew[1] < sublo[1] || xnew[1] > subhi[1]) {
      if (closer_subbox_edge(1, x))
        xnew[1] = subhi[1];
      else
        xnew[1] = sublo[1];
    }
  }
  if (idim <= 1) {
    if (xnew[2] < sublo[2] || xnew[2] > subhi[2]) {
      if (closer_subbox_edge(2, x))
        xnew[2] = subhi[2];
      else
        xnew[2] = sublo[2];
    }
  }
  int proc = point_drop_tiled_recurse(xnew, 0, nprocs - 1);
  if (proc == me) return me;
  if (idim == 0) {
    int done = 1;
    if (rcbinfo[proc].mysplit[1][0] == rcbinfo[me].mysplit[1][1]) {
      xnew[1] -= EPSILON * (subhi[1] - sublo[1]);
      done = 0;
    }
    if (rcbinfo[proc].mysplit[2][0] == rcbinfo[me].mysplit[2][1]) {
      xnew[2] -= EPSILON * (subhi[2] - sublo[2]);
      done = 0;
    }
    if (!done) {
      proc = point_drop_tiled_recurse(xnew, 0, nprocs - 1);
      done = 1;
      if (rcbinfo[proc].mysplit[1][0] == rcbinfo[me].mysplit[1][1]) {
        xnew[1] -= EPSILON * (subhi[1] - sublo[1]);
        done = 0;
      }
      if (rcbinfo[proc].mysplit[2][0] == rcbinfo[me].mysplit[2][1]) {
        xnew[2] -= EPSILON * (subhi[2] - sublo[2]);
        done = 0;
      }
      if (!done) proc = point_drop_tiled_recurse(xnew, 0, nprocs - 1);
    }
  } else if (idim == 1) {
    if (rcbinfo[proc].mysplit[2][0] == rcbinfo[me].mysplit[2][1]) {
      xnew[2] -= EPSILON * (subhi[2] - sublo[2]);
      proc = point_drop_tiled_recurse(xnew, 0, nprocs - 1);
    }
  }
  return proc;
}
int CommTiled::point_drop_tiled_recurse(double *x, int proclower, int procupper)
{
  if (proclower == procupper) return proclower;
  int procmid = proclower + (procupper - proclower) / 2 + 1;
  int idim = rcbinfo[procmid].dim;
  double cut = boxlo[idim] + prd[idim] * rcbinfo[procmid].cutfrac;
  if (x[idim] < cut)
    return point_drop_tiled_recurse(x, proclower, procmid - 1);
  else
    return point_drop_tiled_recurse(x, procmid, procupper);
}
int CommTiled::closer_subbox_edge(int idim, double *x)
{
  double deltalo, deltahi;
  if (sublo[idim] == boxlo[idim])
    deltalo = fabs(x[idim] - prd[idim] - sublo[idim]);
  else
    deltalo = fabs(x[idim] - sublo[idim]);
  if (subhi[idim] == boxhi[idim])
    deltahi = fabs(x[idim] + prd[idim] - subhi[idim]);
  else
    deltahi = fabs(x[idim] - subhi[idim]);
  if (deltalo < deltahi) return 0;
  return 1;
}
void CommTiled::coord2proc_setup()
{
  if (!rcbnew) return;
  if (!rcbinfo) rcbinfo = (RCBinfo *) memory->smalloc(nprocs * sizeof(RCBinfo), "comm:rcbinfo");
  rcbnew = 0;
  RCBinfo rcbone;
  memcpy(&rcbone.mysplit[0][0], &mysplit[0][0], 6 * sizeof(double));
  rcbone.cutfrac = rcbcutfrac;
  rcbone.dim = rcbcutdim;
  MPI_Allgather(&rcbone, sizeof(RCBinfo), MPI_CHAR, rcbinfo, sizeof(RCBinfo), MPI_CHAR, world);
}
int CommTiled::coord2proc(double *x, int &igx, int &igy, int &igz)
{
  if (layout != Comm::LAYOUT_TILED) return Comm::coord2proc(x, igx, igy, igz);
  return point_drop_tiled_recurse(x, 0, nprocs - 1);
}
void CommTiled::grow_send(int n, int flag)
{
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
void CommTiled::grow_recv(int n)
{
  maxrecv = static_cast<int>(BUFFACTOR * n);
  memory->destroy(buf_recv);
  memory->create(buf_recv, maxrecv, "comm:buf_recv");
}
void CommTiled::grow_list(int iswap, int iwhich, int n)
{
  maxsendlist[iswap][iwhich] = static_cast<int>(BUFFACTOR * n);
  memory->grow(sendlist[iswap][iwhich], maxsendlist[iswap][iwhich], "comm:sendlist[i][j]");
}
void CommTiled::allocate_swap(int n)
{
  nsendproc = new int[n];
  nrecvproc = new int[n];
  sendother = new int[n];
  recvother = new int[n];
  sendself = new int[n];
  nprocmax = new int[n];
  sendproc = new int *[n];
  recvproc = new int *[n];
  sendnum = new int *[n];
  recvnum = new int *[n];
  size_forward_recv = new int *[n];
  firstrecv = new int *[n];
  size_reverse_send = new int *[n];
  size_reverse_recv = new int *[n];
  forward_recv_offset = new int *[n];
  reverse_recv_offset = new int *[n];
  pbc_flag = new int *[n];
  pbc = new int **[n];
  sendbox = new double **[n];
  sendbox_multi = new double ***[n];
  sendbox_multiold = new double ***[n];
  maxsendlist = new int *[n];
  sendlist = new int **[n];
  for (int i = 0; i < n; i++) {
    sendproc[i] = recvproc[i] = nullptr;
    sendnum[i] = recvnum[i] = nullptr;
    size_forward_recv[i] = firstrecv[i] = nullptr;
    size_reverse_send[i] = size_reverse_recv[i] = nullptr;
    forward_recv_offset[i] = reverse_recv_offset[i] = nullptr;
    pbc_flag[i] = nullptr;
    pbc[i] = nullptr;
    sendbox[i] = nullptr;
    sendbox_multi[i] = nullptr;
    sendbox_multiold[i] = nullptr;
    maxsendlist[i] = nullptr;
    sendlist[i] = nullptr;
  }
  maxrequest = 0;
  requests = nullptr;
  for (int i = 0; i < n; i++) {
    nprocmax[i] = DELTA_PROCS;
    grow_swap_send(i, DELTA_PROCS, 0);
    grow_swap_recv(i, DELTA_PROCS);
  }
  nexchproc = new int[n / 2];
  nexchprocmax = new int[n / 2];
  exchproc = new int *[n / 2];
  exchnum = new int *[n / 2];
  for (int i = 0; i < n / 2; i++) {
    nexchprocmax[i] = DELTA_PROCS;
    exchproc[i] = new int[DELTA_PROCS];
    exchnum[i] = new int[DELTA_PROCS];
  }
}
void CommTiled::grow_swap_send(int i, int n, int nold)
{
  delete[] sendproc[i];
  sendproc[i] = new int[n];
  delete[] sendnum[i];
  sendnum[i] = new int[n];
  delete[] size_reverse_recv[i];
  size_reverse_recv[i] = new int[n];
  delete[] reverse_recv_offset[i];
  reverse_recv_offset[i] = new int[n];
  delete[] pbc_flag[i];
  pbc_flag[i] = new int[n];
  memory->destroy(pbc[i]);
  memory->create(pbc[i], n, 6, "comm:pbc_flag");
  memory->destroy(sendbox[i]);
  memory->create(sendbox[i], n, 6, "comm:sendbox");
  grow_swap_send_multi(i, n);
  memory->destroy(sendbox_multiold[i]);
  memory->create(sendbox_multiold[i], n, atom->ntypes + 1, 6, "comm:sendbox_multiold");
  delete[] maxsendlist[i];
  maxsendlist[i] = new int[n];
  for (int j = 0; j < nold; j++) memory->destroy(sendlist[i][j]);
  delete[] sendlist[i];
  sendlist[i] = new int *[n];
  for (int j = 0; j < n; j++) {
    maxsendlist[i][j] = BUFMIN;
    memory->create(sendlist[i][j], BUFMIN, "comm:sendlist[i][j]");
  }
}
void CommTiled::grow_swap_recv(int i, int n)
{
  delete[] recvproc[i];
  recvproc[i] = new int[n];
  delete[] recvnum[i];
  recvnum[i] = new int[n];
  delete[] size_forward_recv[i];
  size_forward_recv[i] = new int[n];
  delete[] firstrecv[i];
  firstrecv[i] = new int[n];
  delete[] forward_recv_offset[i];
  forward_recv_offset[i] = new int[n];
  delete[] size_reverse_send[i];
  size_reverse_send[i] = new int[n];
}
void CommTiled::grow_swap_send_multi(int i, int n)
{
  memory->destroy(sendbox_multi[i]);
  if (ncollections > 0) memory->create(sendbox_multi[i], n, ncollections, 6, "comm:sendbox_multi");
}
void CommTiled::deallocate_swap(int n)
{
  delete[] nsendproc;
  delete[] nrecvproc;
  delete[] sendother;
  delete[] recvother;
  delete[] sendself;
  for (int i = 0; i < n; i++) {
    delete[] sendproc[i];
    delete[] recvproc[i];
    delete[] sendnum[i];
    delete[] recvnum[i];
    delete[] size_forward_recv[i];
    delete[] firstrecv[i];
    delete[] size_reverse_send[i];
    delete[] size_reverse_recv[i];
    delete[] forward_recv_offset[i];
    delete[] reverse_recv_offset[i];
    delete[] pbc_flag[i];
    memory->destroy(pbc[i]);
    memory->destroy(sendbox[i]);
    memory->destroy(sendbox_multi[i]);
    memory->destroy(sendbox_multiold[i]);
    delete[] maxsendlist[i];
    for (int j = 0; j < nprocmax[i]; j++) memory->destroy(sendlist[i][j]);
    delete[] sendlist[i];
  }
  delete[] sendproc;
  delete[] recvproc;
  delete[] sendnum;
  delete[] recvnum;
  delete[] size_forward_recv;
  delete[] firstrecv;
  delete[] size_reverse_send;
  delete[] size_reverse_recv;
  delete[] forward_recv_offset;
  delete[] reverse_recv_offset;
  delete[] pbc_flag;
  delete[] pbc;
  delete[] sendbox;
  delete[] sendbox_multi;
  delete[] sendbox_multiold;
  delete[] maxsendlist;
  delete[] sendlist;
  delete[] requests;
  delete[] nprocmax;
  delete[] nexchproc;
  delete[] nexchprocmax;
  for (int i = 0; i < n / 2; i++) {
    delete[] exchproc[i];
    delete[] exchnum[i];
  }
  delete[] exchproc;
  delete[] exchnum;
}
double CommTiled::memory_usage()
{
  double bytes = 0;
  return bytes;
}
