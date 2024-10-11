#include "rcb.h"
#include "memory.h"
#include <cstring>
using namespace LAMMPS_NS;
#define MYHUGE 1.0e30
#define TINY 1.0e-6
#define DELTA 16384
void box_merge(void *, void *, int *, MPI_Datatype *);
void median_merge(void *, void *, int *, MPI_Datatype *);
RCB::RCB(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  ndot = maxdot = 0;
  dots = nullptr;
  nlist = maxlist = 0;
  dotlist = dotmark = dotmark_select = nullptr;
  maxbuf = 0;
  buf = nullptr;
  maxrecv = maxsend = 0;
  recvproc = recvindex = sendproc = sendindex = nullptr;
  tree = nullptr;
  MPI_Type_contiguous(6,MPI_DOUBLE,&box_type);
  MPI_Type_commit(&box_type);
  MPI_Type_contiguous(sizeof(Median),MPI_CHAR,&med_type);
  MPI_Type_commit(&med_type);
  MPI_Op_create(box_merge,1,&box_op);
  MPI_Op_create(median_merge,1,&med_op);
  reuse = 0;
}
RCB::~RCB()
{
  memory->sfree(dots);
  memory->destroy(dotlist);
  memory->destroy(dotmark);
  memory->destroy(dotmark_select);
  memory->sfree(buf);
  memory->destroy(recvproc);
  memory->destroy(recvindex);
  memory->destroy(sendproc);
  memory->destroy(sendindex);
  memory->sfree(tree);
  MPI_Type_free(&med_type);
  MPI_Type_free(&box_type);
  MPI_Op_free(&box_op);
  MPI_Op_free(&med_op);
}
void RCB::compute(int dimension, int n, double **x, double *wt,
                  double *bboxlo, double *bboxhi)
{
  int i,j,k;
  int keep,outgoing,incoming,incoming2;
  int dim,markactive;
  int indexlo,indexhi;
  int first_iteration,breakflag;
  double wttot,wtlo,wthi,wtsum,wtok,wtupto,wtmax;
  double targetlo,targethi;
  double valuemin,valuemax,valuehalf,valuehalf_select,smaller;
  double tolerance;
  MPI_Comm comm,comm_half;
  MPI_Request request,request2;
  Median med,medme;
  ndot = nkeep = noriginal = n;
  if (ndot > maxdot) {
    maxdot = ndot;
    memory->sfree(dots);
    dots = (Dot *) memory->smalloc(ndot*sizeof(Dot),"RCB:dots");
  }
  for (i = 0; i < ndot; i++) {
    dots[i].x[0] = x[i][0];
    dots[i].x[1] = x[i][1];
    dots[i].x[2] = x[i][2];
    dots[i].proc = me;
    dots[i].index = i;
  }
  if (wt)
    for (i = 0; i < ndot; i++) dots[i].wt = wt[i];
  else
    for (i = 0; i < ndot; i++) dots[i].wt = 1.0;
  lo = bbox.lo;
  hi = bbox.hi;
  lo[0] = bboxlo[0];
  lo[1] = bboxlo[1];
  lo[2] = bboxlo[2];
  hi[0] = bboxhi[0];
  hi[1] = bboxhi[1];
  hi[2] = bboxhi[2];
  cut = 0.0;
  cutdim = -1;
  counters[0] = 0;
  counters[1] = 0;
  counters[2] = 0;
  counters[3] = ndot;
  counters[4] = maxdot;
  counters[5] = 0;
  counters[6] = 0;
  MPI_Comm_dup(world,&comm);
  int procpartner,procpartner2;
  int procmid;
  int proclower = 0;
  int procupper = nprocs - 1;
  while (proclower != procupper) {
    procmid = proclower + (procupper - proclower) / 2 + 1;
    if (me < procmid)
      procpartner = me + (procmid - proclower);
    else
      procpartner = me - (procmid - proclower);
    int readnumber = 1;
    if (procpartner > procupper) {
      readnumber = 0;
      procpartner--;
    }
    if (me == procupper && procpartner != procmid - 1) {
      readnumber = 2;
      procpartner2 = procpartner + 1;
    }
    wtmax = wtsum = 0.0;
    if (wt) {
      for (i = 0; i < ndot; i++) {
        wtsum += dots[i].wt;
        if (dots[i].wt > wtmax) wtmax = dots[i].wt;
      }
    } else {
      for (i = 0; i < ndot; i++) wtsum += dots[i].wt;
      wtmax = 1.0;
    }
    MPI_Allreduce(&wtsum,&wttot,1,MPI_DOUBLE,MPI_SUM,comm);
    if (wt) MPI_Allreduce(&wtmax,&tolerance,1,MPI_DOUBLE,MPI_MAX,comm);
    else tolerance = 1.0;
    tolerance *= 1.0 + TINY;
    targetlo = wttot * (procmid - proclower) / (procupper + 1 - proclower);
    targethi = wttot - targetlo;
    int dim_select = -1;
    double largest = -1.0;
    for (dim = 0; dim < dimension; dim++) {
      if (ndot > maxlist) {
        memory->destroy(dotlist);
        memory->destroy(dotmark);
        memory->destroy(dotmark_select);
        maxlist = maxdot;
        memory->create(dotlist,maxlist,"RCB:dotlist");
        memory->create(dotmark,maxlist,"RCB:dotmark");
        memory->create(dotmark_select,maxlist,"RCB:dotmark_select");
      }
      nlist = ndot;
      for (i = 0; i < nlist; i++) dotlist[i] = i;
      wtlo = wthi = 0.0;
      valuemin = lo[dim];
      valuemax = hi[dim];
      first_iteration = 1;
      indexlo = indexhi = 0;
      while (true) {
        if (first_iteration && reuse && dim == tree[procmid].dim) {
          counters[5]++;
          valuehalf = tree[procmid].cut;
          if (valuehalf < valuemin || valuehalf > valuemax)
            valuehalf = 0.5 * (valuemin + valuemax);
        } else if (wt)
          valuehalf = valuemin + (targetlo - wtlo) /
            (wttot - wtlo - wthi) * (valuemax - valuemin);
        else
          valuehalf = 0.5 * (valuemin + valuemax);
        first_iteration = 0;
        medme.totallo = medme.totalhi = 0.0;
        medme.valuelo = -MYHUGE;
        medme.valuehi = MYHUGE;
        medme.wtlo = medme.wthi = 0.0;
        medme.countlo = medme.counthi = 0;
        medme.proclo = medme.prochi = me;
        for (j = 0; j < nlist; j++) {
          i = dotlist[j];
          if (dots[i].x[dim] <= valuehalf) {
            medme.totallo += dots[i].wt;
            dotmark[i] = 0;
            if (dots[i].x[dim] > medme.valuelo) {
              medme.valuelo = dots[i].x[dim];
              medme.wtlo = dots[i].wt;
              medme.countlo = 1;
              indexlo = i;
            } else if (dots[i].x[dim] == medme.valuelo) {
              medme.wtlo += dots[i].wt;
              medme.countlo++;
            }
          }
          else {
            medme.totalhi += dots[i].wt;
            dotmark[i] = 1;
            if (dots[i].x[dim] < medme.valuehi) {
              medme.valuehi = dots[i].x[dim];
              medme.wthi = dots[i].wt;
              medme.counthi = 1;
              indexhi = i;
            } else if (dots[i].x[dim] == medme.valuehi) {
              medme.wthi += dots[i].wt;
              medme.counthi++;
            }
          }
        }
        counters[0]++;
        MPI_Allreduce(&medme,&med,1,med_type,med_op,comm);
        if (wtlo + med.totallo < targetlo) {
          wtlo += med.totallo;
          valuehalf = med.valuehi;
          if (med.counthi == 1) {
            if (wtlo + med.wthi < targetlo) {
              if (me == med.prochi) dotmark[indexhi] = 0;
            }
            else {
              if (wtlo + med.wthi - targetlo < targetlo - wtlo)
                if (me == med.prochi) dotmark[indexhi] = 0;
              break;
            }
          }
          else {
            breakflag = 0;
            wtok = 0.0;
            if (medme.valuehi == med.valuehi) wtok = medme.wthi;
            if (wtlo + med.wthi >= targetlo) {
              MPI_Scan(&wtok,&wtupto,1,MPI_DOUBLE,MPI_SUM,comm);
              wtmax = targetlo - wtlo;
              if (wtupto > wtmax) wtok = wtok - (wtupto - wtmax);
              breakflag = 1;
            }
            for (j = 0, wtsum = 0.0; j < nlist && wtsum < wtok; j++) {
              i = dotlist[j];
              if (dots[i].x[dim] == med.valuehi) {
                if (wtsum + dots[i].wt - wtok < wtok - wtsum)
                  dotmark[i] = 0;
                wtsum += dots[i].wt;
              }
            }
            if (breakflag) break;
          }
          wtlo += med.wthi;
          if (targetlo-wtlo <= tolerance) break;
          valuemin = med.valuehi;
          markactive = 1;
        }
        else if (wthi + med.totalhi < targethi) {
          wthi += med.totalhi;
          valuehalf = med.valuelo;
          if (med.countlo == 1) {
            if (wthi + med.wtlo < targethi) {
              if (me == med.proclo) dotmark[indexlo] = 1;
            }
            else {
              if (wthi + med.wtlo - targethi < targethi - wthi)
                if (me == med.proclo) dotmark[indexlo] = 1;
              break;
            }
          }
          else {
            breakflag = 0;
            wtok = 0.0;
            if (medme.valuelo == med.valuelo) wtok = medme.wtlo;
            if (wthi + med.wtlo >= targethi) {
              MPI_Scan(&wtok,&wtupto,1,MPI_DOUBLE,MPI_SUM,comm);
              wtmax = targethi - wthi;
              if (wtupto > wtmax) wtok = wtok - (wtupto - wtmax);
              breakflag = 1;
            }
            for (j = 0, wtsum = 0.0; j < nlist && wtsum < wtok; j++) {
              i = dotlist[j];
              if (dots[i].x[dim] == med.valuelo) {
                if (wtsum + dots[i].wt - wtok < wtok - wtsum)
                  dotmark[i] = 1;
                wtsum += dots[i].wt;
              }
            }
            if (breakflag) break;
          }
          wthi += med.wtlo;
          if (targethi-wthi <= tolerance) break;
          valuemax = med.valuelo;
          markactive = 0;
        }
        else
          break;
        k = 0;
        for (j = 0; j < nlist; j++) {
          i = dotlist[j];
          if (dotmark[i] == markactive) dotlist[k++] = i;
        }
        nlist = k;
      }
      smaller = MIN(valuehalf-lo[dim],hi[dim]-valuehalf);
      if (smaller > largest) {
        largest = smaller;
        dim_select = dim;
        valuehalf_select = valuehalf;
        if (ndot > 0) memcpy(dotmark_select,dotmark,ndot*sizeof(int));
      }
    }
    dim = dim_select;
    valuehalf = valuehalf_select;
    if (ndot > 0) memcpy(dotmark,dotmark_select,ndot*sizeof(int));
    if (largest == 0.0) {
      dim = 0;
      if (hi[1]-lo[1] > hi[0]-lo[0]) dim = 1;
      if (dimension == 3 && hi[2]-lo[2] > hi[dim]-lo[dim]) dim = 2;
      valuehalf = 0.5* (lo[dim] + hi[dim]);
      for (j = 0; j < nlist; j++) {
        i = dotlist[j];
        if (dots[i].x[dim] <= valuehalf) dotmark[i] = 0;
        else dotmark[i] = 1;
      }
    }
    if (me == procmid) {
      cut = valuehalf;
      cutdim = dim;
    }
    if (me < procmid) hi[dim] = valuehalf;
    else lo[dim] = valuehalf;
    markactive = (me < procpartner);
    for (i = 0, keep = 0, outgoing = 0; i < ndot; i++)
      if (dotmark[i] == markactive) outgoing++;
      else if (i < nkeep) keep++;
    nkeep = keep;
    MPI_Send(&outgoing,1,MPI_INT,procpartner,0,world);
    incoming = 0;
    if (readnumber) {
      MPI_Recv(&incoming,1,MPI_INT,procpartner,0,world,MPI_STATUS_IGNORE);
      if (readnumber == 2) {
        MPI_Recv(&incoming2,1,MPI_INT,procpartner2,0,world,MPI_STATUS_IGNORE);
        incoming += incoming2;
      }
    }
    int ndotnew = ndot - outgoing + incoming;
    if (ndotnew > maxdot) {
      while (maxdot < ndotnew) maxdot += DELTA;
      dots = (Dot *) memory->srealloc(dots,maxdot*sizeof(Dot),"RCB::dots");
      counters[6]++;
    }
    counters[1] += outgoing;
    counters[2] += incoming;
    if (ndotnew > counters[3]) counters[3] = ndotnew;
    if (maxdot > counters[4]) counters[4] = maxdot;
    if (outgoing > maxbuf) {
      memory->sfree(buf);
      maxbuf = outgoing;
      buf = (Dot *) memory->smalloc(maxbuf*sizeof(Dot),"RCB:buf");
    }
    keep = outgoing = 0;
    for (i = 0; i < ndot; i++) {
      if (dotmark[i] == markactive)
        memcpy(&buf[outgoing++],&dots[i],sizeof(Dot));
      else
        memmove(&dots[keep++],&dots[i],sizeof(Dot));
    }
    if (readnumber > 0) {
      MPI_Irecv(&dots[keep],incoming*sizeof(Dot),MPI_CHAR,
                procpartner,1,world,&request);
      if (readnumber == 2) {
        keep += incoming - incoming2;
        MPI_Irecv(&dots[keep],incoming2*sizeof(Dot),MPI_CHAR,
                  procpartner2,1,world,&request2);
      }
    }
    if (readnumber > 0) {
      MPI_Send(nullptr,0,MPI_INT,procpartner,0,world);
      if (readnumber == 2) MPI_Send(nullptr,0,MPI_INT,procpartner2,0,world);
    }
    MPI_Recv(nullptr,0,MPI_INT,procpartner,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(buf,outgoing*sizeof(Dot),MPI_CHAR,procpartner,1,world);
    if (readnumber > 0) {
      MPI_Wait(&request,MPI_STATUS_IGNORE);
      if (readnumber == 2) MPI_Wait(&request2,MPI_STATUS_IGNORE);
    }
    ndot = ndotnew;
    int split;
    if (me < procmid) {
      procupper = procmid - 1;
      split = 0;
    } else {
      proclower = procmid;
      split = 1;
    }
    MPI_Comm_split(comm,split,me,&comm_half);
    MPI_Comm_free(&comm);
    comm = comm_half;
  }
  MPI_Comm_free(&comm);
  nfinal = ndot;
  if (nfinal > maxrecv) {
    memory->destroy(recvproc);
    memory->destroy(recvindex);
    maxrecv = nfinal;
    memory->create(recvproc,maxrecv,"RCB:recvproc");
    memory->create(recvindex,maxrecv,"RCB:recvindex");
  }
  for (i = 0; i < nfinal; i++) {
    recvproc[i] = dots[i].proc;
    recvindex[i] = dots[i].index;
  }
}
void RCB::compute_old(int dimension, int n, double **x, double *wt,
                      double *bboxlo, double *bboxhi)
{
  int i,j,k;
  int keep,outgoing,incoming,incoming2;
  int dim,markactive;
  int indexlo,indexhi;
  int first_iteration,breakflag;
  double wttot,wtlo,wthi,wtsum,wtok,wtupto,wtmax;
  double targetlo,targethi;
  double valuemin,valuemax,valuehalf;
  double tolerance;
  MPI_Comm comm,comm_half;
  MPI_Request request,request2;
  Median med,medme;
  ndot = nkeep = noriginal = n;
  if (ndot > maxdot) {
    maxdot = ndot;
    memory->sfree(dots);
    dots = (Dot *) memory->smalloc(ndot*sizeof(Dot),"RCB:dots");
  }
  for (i = 0; i < ndot; i++) {
    dots[i].x[0] = x[i][0];
    dots[i].x[1] = x[i][1];
    dots[i].x[2] = x[i][2];
    dots[i].proc = me;
    dots[i].index = i;
  }
  if (wt)
    for (i = 0; i < ndot; i++) dots[i].wt = wt[i];
  else
    for (i = 0; i < ndot; i++) dots[i].wt = 1.0;
  lo = bbox.lo;
  hi = bbox.hi;
  lo[0] = bboxlo[0];
  lo[1] = bboxlo[1];
  lo[2] = bboxlo[2];
  hi[0] = bboxhi[0];
  hi[1] = bboxhi[1];
  hi[2] = bboxhi[2];
  cut = 0.0;
  cutdim = -1;
  counters[0] = 0;
  counters[1] = 0;
  counters[2] = 0;
  counters[3] = ndot;
  counters[4] = maxdot;
  counters[5] = 0;
  counters[6] = 0;
  MPI_Comm_dup(world,&comm);
  int procpartner,procpartner2;
  int procmid;
  int proclower = 0;
  int procupper = nprocs - 1;
  while (proclower != procupper) {
    procmid = proclower + (procupper - proclower) / 2 + 1;
    if (me < procmid)
      procpartner = me + (procmid - proclower);
    else
      procpartner = me - (procmid - proclower);
    int readnumber = 1;
    if (procpartner > procupper) {
      readnumber = 0;
      procpartner--;
    }
    if (me == procupper && procpartner != procmid - 1) {
      readnumber = 2;
      procpartner2 = procpartner + 1;
    }
    wtmax = wtsum = 0.0;
    if (wt) {
      for (i = 0; i < ndot; i++) {
        wtsum += dots[i].wt;
        if (dots[i].wt > wtmax) wtmax = dots[i].wt;
      }
    } else {
      for (i = 0; i < ndot; i++) wtsum += dots[i].wt;
      wtmax = 1.0;
    }
    MPI_Allreduce(&wtsum,&wttot,1,MPI_DOUBLE,MPI_SUM,comm);
    if (wt) MPI_Allreduce(&wtmax,&tolerance,1,MPI_DOUBLE,MPI_MAX,comm);
    else tolerance = 1.0;
    tolerance *= 1.0 + TINY;
    targetlo = wttot * (procmid - proclower) / (procupper + 1 - proclower);
    targethi = wttot - targetlo;
    dim = 0;
    if (hi[1]-lo[1] > hi[0]-lo[0]) dim = 1;
    if (dimension == 3) {
      if (dim == 0 && hi[2]-lo[2] > hi[0]-lo[0]) dim = 2;
      if (dim == 1 && hi[2]-lo[2] > hi[1]-lo[1]) dim = 2;
    }
    if (ndot > maxlist) {
      memory->destroy(dotlist);
      memory->destroy(dotmark);
      maxlist = maxdot;
      memory->create(dotlist,maxlist,"RCB:dotlist");
      memory->create(dotmark,maxlist,"RCB:dotmark");
    }
    nlist = ndot;
    for (i = 0; i < nlist; i++) dotlist[i] = i;
    wtlo = wthi = 0.0;
    valuemin = lo[dim];
    valuemax = hi[dim];
    first_iteration = 1;
    indexlo = indexhi = 0;
    while (true) {
      if (first_iteration && reuse && dim == tree[procmid].dim) {
        counters[5]++;
        valuehalf = tree[procmid].cut;
        if (valuehalf < valuemin || valuehalf > valuemax)
          valuehalf = 0.5 * (valuemin + valuemax);
      } else if (wt)
        valuehalf = valuemin + (targetlo - wtlo) /
          (wttot - wtlo - wthi) * (valuemax - valuemin);
      else
        valuehalf = 0.5 * (valuemin + valuemax);
      first_iteration = 0;
      medme.totallo = medme.totalhi = 0.0;
      medme.valuelo = -MYHUGE;
      medme.valuehi = MYHUGE;
      medme.wtlo = medme.wthi = 0.0;
      medme.countlo = medme.counthi = 0;
      medme.proclo = medme.prochi = me;
      for (j = 0; j < nlist; j++) {
        i = dotlist[j];
        if (dots[i].x[dim] <= valuehalf) {
          medme.totallo += dots[i].wt;
          dotmark[i] = 0;
          if (dots[i].x[dim] > medme.valuelo) {
            medme.valuelo = dots[i].x[dim];
            medme.wtlo = dots[i].wt;
            medme.countlo = 1;
            indexlo = i;
          } else if (dots[i].x[dim] == medme.valuelo) {
            medme.wtlo += dots[i].wt;
            medme.countlo++;
          }
        }
        else {
          medme.totalhi += dots[i].wt;
          dotmark[i] = 1;
          if (dots[i].x[dim] < medme.valuehi) {
            medme.valuehi = dots[i].x[dim];
            medme.wthi = dots[i].wt;
            medme.counthi = 1;
            indexhi = i;
          } else if (dots[i].x[dim] == medme.valuehi) {
            medme.wthi += dots[i].wt;
            medme.counthi++;
          }
        }
      }
      counters[0]++;
      MPI_Allreduce(&medme,&med,1,med_type,med_op,comm);
      if (wtlo + med.totallo < targetlo) {
        wtlo += med.totallo;
        valuehalf = med.valuehi;
        if (med.counthi == 1) {
          if (wtlo + med.wthi < targetlo) {
            if (me == med.prochi) dotmark[indexhi] = 0;
          }
          else {
            if (wtlo + med.wthi - targetlo < targetlo - wtlo)
              if (me == med.prochi) dotmark[indexhi] = 0;
            break;
          }
        }
        else {
          breakflag = 0;
          wtok = 0.0;
          if (medme.valuehi == med.valuehi) wtok = medme.wthi;
          if (wtlo + med.wthi >= targetlo) {
            MPI_Scan(&wtok,&wtupto,1,MPI_DOUBLE,MPI_SUM,comm);
            wtmax = targetlo - wtlo;
            if (wtupto > wtmax) wtok = wtok - (wtupto - wtmax);
            breakflag = 1;
          }
          for (j = 0, wtsum = 0.0; j < nlist && wtsum < wtok; j++) {
            i = dotlist[j];
            if (dots[i].x[dim] == med.valuehi) {
              if (wtsum + dots[i].wt - wtok < wtok - wtsum)
                dotmark[i] = 0;
              wtsum += dots[i].wt;
            }
          }
          if (breakflag) break;
        }
        wtlo += med.wthi;
        if (targetlo-wtlo <= tolerance) break;
        valuemin = med.valuehi;
        markactive = 1;
      }
      else if (wthi + med.totalhi < targethi) {
        wthi += med.totalhi;
        valuehalf = med.valuelo;
        if (med.countlo == 1) {
          if (wthi + med.wtlo < targethi) {
            if (me == med.proclo) dotmark[indexlo] = 1;
          }
          else {
            if (wthi + med.wtlo - targethi < targethi - wthi)
              if (me == med.proclo) dotmark[indexlo] = 1;
            break;
          }
        }
        else {
          breakflag = 0;
          wtok = 0.0;
          if (medme.valuelo == med.valuelo) wtok = medme.wtlo;
          if (wthi + med.wtlo >= targethi) {
            MPI_Scan(&wtok,&wtupto,1,MPI_DOUBLE,MPI_SUM,comm);
            wtmax = targethi - wthi;
            if (wtupto > wtmax) wtok = wtok - (wtupto - wtmax);
            breakflag = 1;
          }
          for (j = 0, wtsum = 0.0; j < nlist && wtsum < wtok; j++) {
            i = dotlist[j];
            if (dots[i].x[dim] == med.valuelo) {
              if (wtsum + dots[i].wt - wtok < wtok - wtsum)
                dotmark[i] = 1;
              wtsum += dots[i].wt;
            }
          }
          if (breakflag) break;
        }
        wthi += med.wtlo;
        if (targethi-wthi <= tolerance) break;
        valuemax = med.valuelo;
        markactive = 0;
      }
      else
        break;
      k = 0;
      for (j = 0; j < nlist; j++) {
        i = dotlist[j];
        if (dotmark[i] == markactive) dotlist[k++] = i;
      }
      nlist = k;
    }
    if (me == procmid) {
      cut = valuehalf;
      cutdim = dim;
    }
    if (me < procmid) hi[dim] = valuehalf;
    else lo[dim] = valuehalf;
    markactive = (me < procpartner);
    for (i = 0, keep = 0, outgoing = 0; i < ndot; i++)
      if (dotmark[i] == markactive) outgoing++;
      else if (i < nkeep) keep++;
    nkeep = keep;
    MPI_Send(&outgoing,1,MPI_INT,procpartner,0,world);
    incoming = 0;
    if (readnumber) {
      MPI_Recv(&incoming,1,MPI_INT,procpartner,0,world,MPI_STATUS_IGNORE);
      if (readnumber == 2) {
        MPI_Recv(&incoming2,1,MPI_INT,procpartner2,0,world,MPI_STATUS_IGNORE);
        incoming += incoming2;
      }
    }
    int ndotnew = ndot - outgoing + incoming;
    if (ndotnew > maxdot) {
      while (maxdot < ndotnew) maxdot += DELTA;
      dots = (Dot *) memory->srealloc(dots,maxdot*sizeof(Dot),"RCB::dots");
      counters[6]++;
    }
    counters[1] += outgoing;
    counters[2] += incoming;
    if (ndotnew > counters[3]) counters[3] = ndotnew;
    if (maxdot > counters[4]) counters[4] = maxdot;
    if (outgoing > maxbuf) {
      memory->sfree(buf);
      maxbuf = outgoing;
      buf = (Dot *) memory->smalloc(maxbuf*sizeof(Dot),"RCB:buf");
    }
    keep = outgoing = 0;
    for (i = 0; i < ndot; i++) {
      if (dotmark[i] == markactive)
        memcpy(&buf[outgoing++],&dots[i],sizeof(Dot));
      else
        memmove(&dots[keep++],&dots[i],sizeof(Dot));
    }
    if (readnumber > 0) {
      MPI_Irecv(&dots[keep],incoming*sizeof(Dot),MPI_CHAR,
                procpartner,1,world,&request);
      if (readnumber == 2) {
        keep += incoming - incoming2;
        MPI_Irecv(&dots[keep],incoming2*sizeof(Dot),MPI_CHAR,
                  procpartner2,1,world,&request2);
      }
    }
    if (readnumber > 0) {
      MPI_Send(nullptr,0,MPI_INT,procpartner,0,world);
      if (readnumber == 2) MPI_Send(nullptr,0,MPI_INT,procpartner2,0,world);
    }
    MPI_Recv(nullptr,0,MPI_INT,procpartner,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(buf,outgoing*sizeof(Dot),MPI_CHAR,procpartner,1,world);
    if (readnumber > 0) {
      MPI_Wait(&request,MPI_STATUS_IGNORE);
      if (readnumber == 2) MPI_Wait(&request2,MPI_STATUS_IGNORE);
    }
    ndot = ndotnew;
    int split;
    if (me < procmid) {
      procupper = procmid - 1;
      split = 0;
    } else {
      proclower = procmid;
      split = 1;
    }
    MPI_Comm_split(comm,split,me,&comm_half);
    MPI_Comm_free(&comm);
    comm = comm_half;
  }
  MPI_Comm_free(&comm);
  nfinal = ndot;
  if (nfinal > maxrecv) {
    memory->destroy(recvproc);
    memory->destroy(recvindex);
    maxrecv = nfinal;
    memory->create(recvproc,maxrecv,"RCB:recvproc");
    memory->create(recvindex,maxrecv,"RCB:recvindex");
  }
  for (i = 0; i < nfinal; i++) {
    recvproc[i] = dots[i].proc;
    recvindex[i] = dots[i].index;
  }
}
void box_merge(void *in, void *inout, int * , MPI_Datatype * )
{
  auto box1 = (RCB::BBox *) in;
  auto box2 = (RCB::BBox *) inout;
  for (int i = 0; i < 3; i++) {
    if (box1->lo[i] < box2->lo[i]) box2->lo[i] = box1->lo[i];
    if (box1->hi[i] > box2->hi[i]) box2->hi[i] = box1->hi[i];
  }
}
void median_merge(void *in, void *inout, int * , MPI_Datatype * )
{
  auto med1 = (RCB::Median *) in;
  auto med2 = (RCB::Median *) inout;
  med2->totallo += med1->totallo;
  if (med1->valuelo > med2->valuelo) {
    med2->valuelo = med1->valuelo;
    med2->wtlo = med1->wtlo;
    med2->countlo = med1->countlo;
    med2->proclo = med1->proclo;
  }
  else if (med1->valuelo == med2->valuelo) {
    med2->wtlo += med1->wtlo;
    med2->countlo += med1->countlo;
    if (med1->proclo < med2->proclo) med2->proclo = med1->proclo;
  }
  med2->totalhi += med1->totalhi;
  if (med1->valuehi < med2->valuehi) {
    med2->valuehi = med1->valuehi;
    med2->wthi = med1->wthi;
    med2->counthi = med1->counthi;
    med2->prochi = med1->prochi;
  }
  else if (med1->valuehi == med2->valuehi) {
    med2->wthi += med1->wthi;
    med2->counthi += med1->counthi;
    if (med1->prochi < med2->prochi) med2->prochi = med1->prochi;
  }
}
void RCB::invert(int sortflag)
{
  int nsend = nfinal-nkeep;
  int *proclist;
  memory->create(proclist,nsend,"RCB:proclist");
  auto sinvert = (Invert *) memory->smalloc(nsend*sizeof(Invert),"RCB:sinvert");
  int m = 0;
  for (int i = nkeep; i < nfinal; i++) {
    proclist[m] = recvproc[i];
    sinvert[m].rindex = recvindex[i];
    sinvert[m].sproc = me;
    sinvert[m].sindex = i;
    m++;
  }
  if (noriginal > maxsend) {
    memory->destroy(sendproc);
    memory->destroy(sendindex);
    maxsend = noriginal;
    memory->create(sendproc,maxsend,"RCB:sendproc");
    memory->create(sendindex,maxsend,"RCB:sendindex");
  }
  for (int i = 0; i < nkeep; i++) {
    sendproc[recvindex[i]] = me;
    sendindex[recvindex[i]] = i;
  }
  memory->destroy(proclist);
  memory->destroy(sinvert);
}
double RCB::memory_usage()
{
  double bytes = 0;
  return bytes;
}
