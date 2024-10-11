#include "balance.h"
#include "update.h"
#include "atom.h"
#include "neighbor.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "imbalance.h"
#include "imbalance_group.h"
#include "imbalance_neigh.h"
#include "imbalance_store.h"
#include "imbalance_var.h"
#include "memory.h"
#include "modify.h"
#include "pair.h"
#include "rcb.h"
#include "error.h"
#include <cmath>
#include <cstring>
using namespace LAMMPS_NS;
double EPSNEIGH = 1.0e-3;
enum{XYZ,SHIFT,BISECTION};
enum{NONE,UNIFORM,USER};
enum{X,Y,Z};
Balance::Balance(LAMMPS *lmp) : Command(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  user_xsplit = user_ysplit = user_zsplit = nullptr;
  shift_allocate = 0;
  proccost = allproccost = nullptr;
  rcb = nullptr;
  nimbalance = 0;
  imbalances = nullptr;
  fp = nullptr;
  firststep = 1;
}
Balance::~Balance()
{
  memory->destroy(proccost);
  memory->destroy(allproccost);
  delete[] user_xsplit;
  delete[] user_ysplit;
  delete[] user_zsplit;
  if (shift_allocate) {
    delete[] bdim;
    delete[] onecost;
    delete[] allcost;
    delete[] sum;
    delete[] target;
    delete[] lo;
    delete[] hi;
    delete[] losum;
    delete[] hisum;
  }
  delete rcb;
  for (int i = 0; i < nimbalance; i++) delete imbalances[i];
  delete[] imbalances;
  if (fp) fclose(fp);
}
void Balance::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Balance command before simulation box is defined");
  if (me == 0) utils::logmesg(lmp,"Balancing ...\n");
  if (narg < 2) error->all(FLERR,"Illegal balance command");
  thresh = utils::numeric(FLERR,arg[0],false,lmp);
  int dimension = domain->dimension;
  int *procgrid = comm->procgrid;
  style = -1;
  xflag = yflag = zflag = NONE;
  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"x") == 0) {
      if (style != -1 && style != XYZ)
        error->all(FLERR,"Illegal balance command");
      style = XYZ;
      if (strcmp(arg[iarg+1],"uniform") == 0) {
        if (iarg+2 > narg) error->all(FLERR,"Illegal balance command");
        xflag = UNIFORM;
        iarg += 2;
      } else {
        if (1 + procgrid[0]-1 > narg)
          error->all(FLERR,"Illegal balance command");
        xflag = USER;
        delete[] user_xsplit;
        user_xsplit = new double[procgrid[0]+1];
        user_xsplit[0] = 0.0;
        iarg++;
        for (int i = 1; i < procgrid[0]; i++)
          user_xsplit[i] = utils::numeric(FLERR,arg[iarg++],false,lmp);
        user_xsplit[procgrid[0]] = 1.0;
      }
    } else if (strcmp(arg[iarg],"y") == 0) {
      if (style != -1 && style != XYZ)
        error->all(FLERR,"Illegal balance command");
      style = XYZ;
      if (strcmp(arg[iarg+1],"uniform") == 0) {
        if (iarg+2 > narg) error->all(FLERR,"Illegal balance command");
        yflag = UNIFORM;
        iarg += 2;
      } else {
        if (1 + procgrid[1]-1 > narg)
          error->all(FLERR,"Illegal balance command");
        yflag = USER;
        delete[] user_ysplit;
        user_ysplit = new double[procgrid[1]+1];
        user_ysplit[0] = 0.0;
        iarg++;
        for (int i = 1; i < procgrid[1]; i++)
          user_ysplit[i] = utils::numeric(FLERR,arg[iarg++],false,lmp);
        user_ysplit[procgrid[1]] = 1.0;
      }
    } else if (strcmp(arg[iarg],"z") == 0) {
      if (style != -1 && style != XYZ)
        error->all(FLERR,"Illegal balance command");
      style = XYZ;
      if (strcmp(arg[iarg+1],"uniform") == 0) {
        if (iarg+2 > narg) error->all(FLERR,"Illegal balance command");
        zflag = UNIFORM;
        iarg += 2;
      } else {
        if (1 + procgrid[2]-1 > narg)
          error->all(FLERR,"Illegal balance command");
        zflag = USER;
        delete[] user_zsplit;
        user_zsplit = new double[procgrid[2]+1];
        user_zsplit[0] = 0.0;
        iarg++;
        for (int i = 1; i < procgrid[2]; i++)
          user_zsplit[i] = utils::numeric(FLERR,arg[iarg++],false,lmp);
        user_zsplit[procgrid[2]] = 1.0;
      }
    } else if (strcmp(arg[iarg],"shift") == 0) {
      if (style != -1) error->all(FLERR,"Illegal balance command");
      if (iarg+4 > narg) error->all(FLERR,"Illegal balance command");
      style = SHIFT;
      if (strlen(arg[iarg+1]) > BSTR_SIZE) error->all(FLERR,"Illegal balance command");
      strncpy(bstr,arg[iarg+1],BSTR_SIZE+1);
      nitermax = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
      if (nitermax <= 0) error->all(FLERR,"Illegal balance command");
      stopthresh = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      if (stopthresh < 1.0) error->all(FLERR,"Illegal balance command");
      iarg += 4;
    } else if (strcmp(arg[iarg],"rcb") == 0) {
      if (style != -1) error->all(FLERR,"Illegal balance command");
      style = BISECTION;
      iarg++;
    } else break;
  }
  if (style == XYZ) {
    if (zflag != NONE && dimension == 2)
      error->all(FLERR,"Cannot balance in z dimension for 2d simulation");
    if (xflag == USER)
      for (int i = 1; i <= procgrid[0]; i++)
        if (user_xsplit[i-1] >= user_xsplit[i])
          error->all(FLERR,"Illegal balance command");
    if (yflag == USER)
      for (int i = 1; i <= procgrid[1]; i++)
        if (user_ysplit[i-1] >= user_ysplit[i])
          error->all(FLERR,"Illegal balance command");
    if (zflag == USER)
      for (int i = 1; i <= procgrid[2]; i++)
        if (user_zsplit[i-1] >= user_zsplit[i])
          error->all(FLERR,"Illegal balance command");
  }
  if (style == SHIFT) {
    const int blen=strlen(bstr);
    for (int i = 0; i < blen; i++) {
      if (bstr[i] != 'x' && bstr[i] != 'y' && bstr[i] != 'z')
        error->all(FLERR,"Balance shift string is invalid");
      if (bstr[i] == 'z' && dimension == 2)
        error->all(FLERR,"Balance shift string is invalid");
      for (int j = i+1; j < blen; j++)
        if (bstr[i] == bstr[j])
          error->all(FLERR,"Balance shift string is invalid");
    }
  }
  if (style == BISECTION && comm->style == Comm::BRICK)
    error->all(FLERR,"Balance rcb cannot be used with comm_style brick");
  options(iarg,narg,arg,1);
  MPI_Barrier(world);
  double start_time = platform::walltime();
  lmp->init();
  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  comm->setup();
  comm->exchange();
  if (atom->map_style != Atom::MAP_NONE) atom->map_set();
  if (domain->triclinic) domain->lamda2x(atom->nlocal);
  double maxinit;
  double imbinit = imbalance_factor(maxinit);
  if (comm->layout == Comm::LAYOUT_TILED && style != BISECTION) {
  } else if (imbinit < thresh) return;
#ifdef BALANCE_DEBUG
  if (outflag) dumpout(update->ntimestep);
#endif
  int niter = 0;
  if (style == XYZ) {
    if (comm->layout == Comm::LAYOUT_UNIFORM) {
      if (xflag == USER || yflag == USER || zflag == USER)
        comm->layout = Comm::LAYOUT_NONUNIFORM;
    } else if (comm->layout == Comm::LAYOUT_NONUNIFORM) {
      if (xflag == UNIFORM && yflag == UNIFORM && zflag == UNIFORM)
        comm->layout = Comm::LAYOUT_UNIFORM;
    } else if (comm->layout == Comm::LAYOUT_TILED) {
      if (xflag == UNIFORM && yflag == UNIFORM && zflag == UNIFORM)
        comm->layout = Comm::LAYOUT_UNIFORM;
      else comm->layout = Comm::LAYOUT_NONUNIFORM;
    }
    if (xflag == UNIFORM) {
      for (int i = 0; i < procgrid[0]; i++)
        comm->xsplit[i] = i * 1.0/procgrid[0];
      comm->xsplit[procgrid[0]] = 1.0;
    } else if (xflag == USER)
      for (int i = 0; i <= procgrid[0]; i++) comm->xsplit[i] = user_xsplit[i];
    if (yflag == UNIFORM) {
      for (int i = 0; i < procgrid[1]; i++)
        comm->ysplit[i] = i * 1.0/procgrid[1];
      comm->ysplit[procgrid[1]] = 1.0;
    } else if (yflag == USER)
      for (int i = 0; i <= procgrid[1]; i++) comm->ysplit[i] = user_ysplit[i];
    if (zflag == UNIFORM) {
      for (int i = 0; i < procgrid[2]; i++)
        comm->zsplit[i] = i * 1.0/procgrid[2];
      comm->zsplit[procgrid[2]] = 1.0;
    } else if (zflag == USER)
      for (int i = 0; i <= procgrid[2]; i++) comm->zsplit[i] = user_zsplit[i];
  }
  if (style == SHIFT) {
    comm->layout = Comm::LAYOUT_NONUNIFORM;
    shift_setup_static(bstr);
    niter = shift();
  }
  if (style == BISECTION) {
    comm->layout = Comm::LAYOUT_TILED;
    bisection();
  }
  if (domain->triclinic) domain->set_lamda_box();
  domain->set_local_box();
  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  if (domain->triclinic) domain->lamda2x(atom->nlocal);
  if (outflag) dumpout(update->ntimestep);
  modify->reset_grid();
  if (force->pair) force->pair->reset_grid();
  bigint natoms;
  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (natoms != atom->natoms)
    error->all(FLERR,"Lost atoms via balance: original {}  current {}",
               atom->natoms,natoms);
  double maxfinal;
  double imbfinal = imbalance_factor(maxfinal);
  if (me == 0) {
    std::string mesg = fmt::format(" rebalancing time: {:.3f} seconds\n",
                                   platform::walltime()-start_time);
    mesg += fmt::format("  iteration count = {}\n",niter);
    for (int i = 0; i < nimbalance; ++i) mesg += imbalances[i]->info();
    mesg += fmt::format("  initial/final maximal load/proc = {:.8} {:.8}\n"
                        "  initial/final imbalance factor  = {:.8} {:.8}\n",
                        maxinit,maxfinal,imbinit,imbfinal);
    if (style != BISECTION) {
      mesg += "  x cuts:";
      for (int i = 0; i <= comm->procgrid[0]; i++)
        mesg += fmt::format(" {:.8}",comm->xsplit[i]);
      mesg += "\n  y cuts:";
      for (int i = 0; i <= comm->procgrid[1]; i++)
        mesg += fmt::format(" {:.8}",comm->ysplit[i]);
      mesg += "\n  z cuts:";
      for (int i = 0; i <= comm->procgrid[2]; i++)
        mesg += fmt::format(" {:.8}",comm->zsplit[i]);
      mesg += "\n";
    }
    utils::logmesg(lmp,mesg);
  }
}
void Balance::options(int iarg, int narg, char **arg, int sortflag_default)
{
  nimbalance = 0;
  for (int i = iarg; i < narg; i++)
    if (strcmp(arg[i],"weight") == 0) nimbalance++;
  if (nimbalance) imbalances = new Imbalance*[nimbalance];
  nimbalance = 0;
  sortflag = sortflag_default;
  outflag = 0;
  int outarg = 0;
  fp = nullptr;
  oldrcb = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg+1],"sort") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "balance sort", error);
      sortflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"out") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal (fix) balance command");
      outflag = 1;
      outarg = iarg+1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"old") == 0) {
      oldrcb = 1;
      iarg++;
    } else error->all(FLERR,"Illegal (fix) balance command");
  }
  if (outflag && comm->me == 0) {
    fp = fopen(arg[outarg],"w");
    if (fp == nullptr)
      error->one(FLERR,"Cannot open (fix) balance output file {}: {}",
                                   arg[outarg], utils::getsyserror());
  }
}
double Balance::imbalance_factor(double &maxcost)
{
  double mycost,totalcost;
  mycost = atom->nlocal;
  MPI_Allreduce(&mycost,&maxcost,1,MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(&mycost,&totalcost,1,MPI_DOUBLE,MPI_SUM,world);
  double imbalance = 1.0;
  if (maxcost > 0.0) imbalance = maxcost / (totalcost/nprocs);
  return imbalance;
}
int *Balance::bisection()
{
  if (!rcb) rcb = new RCB(lmp);
  int dim = domain->dimension;
  int triclinic = domain->triclinic;
  double *boxlo,*boxhi,*prd;
  if (triclinic == 0) {
    boxlo = domain->boxlo;
    boxhi = domain->boxhi;
    prd = domain->prd;
  } else {
    boxlo = domain->boxlo_lamda;
    boxhi = domain->boxhi_lamda;
    prd = domain->prd_lamda;
  }
  double shrink[6],shrinkall[6];
  shrink[0] = boxhi[0]; shrink[1] = boxhi[1]; shrink[2] = boxhi[2];
  shrink[3] = boxlo[0]; shrink[4] = boxlo[1]; shrink[5] = boxlo[2];
  double **x = atom->x;
  int nlocal = atom->nlocal;
  if (triclinic) domain->x2lamda(nlocal);
  for (int i = 0; i < nlocal; i++) {
    shrink[0] = MIN(shrink[0],x[i][0]);
    shrink[1] = MIN(shrink[1],x[i][1]);
    shrink[2] = MIN(shrink[2],x[i][2]);
    shrink[3] = MAX(shrink[3],x[i][0]);
    shrink[4] = MAX(shrink[4],x[i][1]);
    shrink[5] = MAX(shrink[5],x[i][2]);
  }
  shrink[3] = -shrink[3]; shrink[4] = -shrink[4]; shrink[5] = -shrink[5];
  MPI_Allreduce(shrink,shrinkall,6,MPI_DOUBLE,MPI_MIN,world);
  shrinkall[3] = -shrinkall[3];
  shrinkall[4] = -shrinkall[4];
  shrinkall[5] = -shrinkall[5];
  double *shrinklo = &shrinkall[0];
  double *shrinkhi = &shrinkall[3];
  if (shrinklo[0] == shrinkhi[0]) {
    shrinklo[0] = boxlo[0];
    shrinkhi[0] = boxhi[0];
  }
  if (shrinklo[1] == shrinkhi[1]) {
    shrinklo[1] = boxlo[1];
    shrinkhi[1] = boxhi[1];
  }
  if (shrinklo[2] == shrinkhi[2]) {
    shrinklo[2] = boxlo[2];
    shrinkhi[2] = boxhi[2];
  }
  if (oldrcb) {
    rcb->compute_old(dim,atom->nlocal,atom->x,nullptr,shrinklo,shrinkhi);
  } else {
    rcb->compute(dim,atom->nlocal,atom->x,nullptr,shrinklo,shrinkhi);
  }
  if (triclinic) domain->lamda2x(nlocal);
  rcb->invert(sortflag);
  double *lo = rcb->lo;
  double *hi = rcb->hi;
  if (lo[0] == shrinklo[0]) lo[0] = boxlo[0];
  if (lo[1] == shrinklo[1]) lo[1] = boxlo[1];
  if (lo[2] == shrinklo[2]) lo[2] = boxlo[2];
  if (hi[0] == shrinkhi[0]) hi[0] = boxhi[0];
  if (hi[1] == shrinkhi[1]) hi[1] = boxhi[1];
  if (hi[2] == shrinkhi[2]) hi[2] = boxhi[2];
  comm->rcbnew = 1;
  int idim = rcb->cutdim;
  if (idim >= 0) comm->rcbcutfrac = (rcb->cut - boxlo[idim]) / prd[idim];
  else comm->rcbcutfrac = 0.0;
  comm->rcbcutdim = idim;
  double (*mysplit)[2] = comm->mysplit;
  mysplit[0][0] = (lo[0] - boxlo[0]) / prd[0];
  if (hi[0] == boxhi[0]) mysplit[0][1] = 1.0;
  else mysplit[0][1] = (hi[0] - boxlo[0]) / prd[0];
  mysplit[1][0] = (lo[1] - boxlo[1]) / prd[1];
  if (hi[1] == boxhi[1]) mysplit[1][1] = 1.0;
  else mysplit[1][1] = (hi[1] - boxlo[1]) / prd[1];
  mysplit[2][0] = (lo[2] - boxlo[2]) / prd[2];
  if (hi[2] == boxhi[2]) mysplit[2][1] = 1.0;
  else mysplit[2][1] = (hi[2] - boxlo[2]) / prd[2];
  return rcb->sendproc;
}
void Balance::shift_setup_static(char *str)
{
  shift_allocate = 1;
  memory->create(proccost,nprocs,"balance:proccost");
  memory->create(allproccost,nprocs,"balance:allproccost");
  ndim = strlen(str);
  bdim = new int[ndim];
  for (int i = 0; i < ndim; i++) {
    if (str[i] == 'x') bdim[i] = X;
    if (str[i] == 'y') bdim[i] = Y;
    if (str[i] == 'z') bdim[i] = Z;
  }
  int max = MAX(comm->procgrid[0],comm->procgrid[1]);
  max = MAX(max,comm->procgrid[2]);
  onecost = new double[max];
  allcost = new double[max];
  sum = new double[max+1];
  target = new double[max+1];
  lo = new double[max+1];
  hi = new double[max+1];
  losum = new double[max+1];
  hisum = new double[max+1];
  if (comm->layout == Comm::LAYOUT_TILED) {
    int *procgrid = comm->procgrid;
    double *xsplit = comm->xsplit;
    double *ysplit = comm->ysplit;
    double *zsplit = comm->zsplit;
    for (int i = 0; i < procgrid[0]; i++) xsplit[i] = i * 1.0/procgrid[0];
    for (int i = 0; i < procgrid[1]; i++) ysplit[i] = i * 1.0/procgrid[1];
    for (int i = 0; i < procgrid[2]; i++) zsplit[i] = i * 1.0/procgrid[2];
    xsplit[procgrid[0]] = ysplit[procgrid[1]] = zsplit[procgrid[2]] = 1.0;
  }
  rho = 0;
}
void Balance::shift_setup(char *str, int nitermax_in, double thresh_in)
{
  shift_setup_static(str);
  nitermax = nitermax_in;
  stopthresh = thresh_in;
  rho = 1;
}
int Balance::shift()
{
  int i,j,k,m,np;
  double mycost,totalcost,boxsize;
  double *split;
  bigint natoms = atom->natoms;
  if (natoms == 0) return 0;
  double delta = pow(stopthresh,1.0/ndim) - 1.0;
  int *procgrid = comm->procgrid;
  domain->x2lamda(atom->nlocal);
  double *prd = domain->prd;
  int niter = 0;
  for (int idim = 0; idim < ndim; idim++) {
    if (bdim[idim] == X) {
      split = comm->xsplit;
      boxsize = prd[0];
    } else if (bdim[idim] == Y) {
      split = comm->ysplit;
      boxsize = prd[1];
    } else if (bdim[idim] == Z) {
      split = comm->zsplit;
      boxsize = prd[2];
    } else continue;
    np = procgrid[bdim[idim]];
    tally(bdim[idim],np,split);
    mycost = atom->nlocal;
    MPI_Allreduce(&mycost,&totalcost,1,MPI_DOUBLE,MPI_SUM,world);
    for (i = 0; i < np; i++) target[i] = totalcost/np * i;
    target[np] = totalcost;
    lo[0] = hi[0] = 0.0;
    lo[np] = hi[np] = 1.0;
    losum[0] = hisum[0] = 0.0;
    losum[np] = hisum[np] = totalcost;
    for (i = 1; i < np; i++) {
      for (j = i; j >= 0; j--)
        if (sum[j] <= target[i]) {
          lo[i] = split[j];
          losum[i] = sum[j];
          break;
        }
      for (j = i; j <= np; j++)
        if (sum[j] >= target[i]) {
          hi[i] = split[j];
          hisum[i] = sum[j];
          break;
        }
    }
#ifdef BALANCE_DEBUG
    if (me == 0) debug_shift_output(idim,0,np,split);
#endif
    int doneflag;
    int change = 1;
    for (m = 0; m < nitermax; m++) {
      change = adjust(np,split);
      tally(bdim[idim],np,split);
      niter++;
#ifdef BALANCE_DEBUG
      if (me == 0) debug_shift_output(idim,m+1,np,split);
      if (outflag) dumpout(update->ntimestep);
#endif
      if (!change) break;
      doneflag = 1;
      for (i = 1; i < np; i++)
        if (fabs(1.0*(sum[i]-target[i]))/target[i] > delta) doneflag = 0;
      if (doneflag) break;
    }
    int duplicate = 0;
    for (i = 1; i < np-1; i++)
      if (split[i] == split[i+1]) duplicate = 1;
    if (duplicate) {
      for (i = 0; i < np; i++)
        lo[i] = 0.5 * (split[i] + split[i+1]);
      i = 1;
      while (i < np-1) {
        j = i+1;
        while (split[j] == split[i]) j++;
        j--;
        if (j > i) {
          delta = (lo[j] - lo[i-1]) / (j-i+2);
          for (k = i; k <= j; k++)
            split[k] = lo[i-1] + (k-i+1)*delta;
        }
        i = j+1;
      }
    }
    double close = (1.0+EPSNEIGH) * neighbor->skin / boxsize;
    double midpt,start,stop,lbound,ubound,spacing;
    i = 0;
    while (i < np) {
      if (split[i+1] - split[i] < close) {
        j = i+1;
        while (true) {
          delta = (j-i) * close;
          midpt = 0.5 * (split[i]+split[j]);
          start = midpt - 0.5*delta;
          stop = midpt + 0.5*delta;
          if (i > 0) lbound = split[i-1] + close;
          else lbound = 0.0;
          if (j < np) ubound = split[j+1] - close;
          else ubound = 1.0;
          if (start >= lbound && stop <= ubound) break;
          if (start < lbound) {
            start = lbound;
            stop = start + delta;
            if (stop <= ubound) break;
          } else if (stop > ubound) {
            stop = ubound;
            start = stop - delta;
            if (start >= lbound) break;
          }
          if (i == 0 && j == np) {
            start = 0.0; stop = 1.0;
            break;
          }
          if (i == 0) j++;
          else if (j == np) i--;
          else if (midpt-lbound < ubound-midpt) i--;
          else j++;
        }
        spacing = (stop-start) / (j-i);
        for (m = i; m <= j; m++)
          split[m] = start + (m-i)*spacing;
        if (j == np) split[np] = 1.0;
        i = j+1;
      } else i++;
    }
    int bad = 0;
    for (i = 0; i < np; i++)
      if (split[i] >= split[i+1]) bad = 1;
    if (bad) error->all(FLERR,"Balance produced bad splits");
    double imbfactor = imbalance_splits();
    if (imbfactor <= stopthresh) break;
  }
  domain->lamda2x(atom->nlocal);
  return niter;
}
void Balance::tally(int dim, int n, double *split)
{
  for (int i = 0; i < n; i++) onecost[i] = 0.0;
  double **x = atom->x;
  int nlocal = atom->nlocal;
  int index;
    for (int i = 0; i < nlocal; i++) {
      index = utils::binary_search(x[i][dim],n,split);
      onecost[index] += 1.0;
    }
  MPI_Allreduce(onecost,allcost,n,MPI_DOUBLE,MPI_SUM,world);
  sum[0] = 0.0;
  for (int i = 1; i < n+1; i++)
    sum[i] = sum[i-1] + allcost[i-1];
}
int Balance::adjust(int n, double *split)
{
  int i;
  double fraction;
  for (i = 1; i < n; i++) {
    if (sum[i] <= target[i]) {
      lo[i] = split[i];
      losum[i] = sum[i];
    }
    if (sum[i] >= target[i]) {
      hi[i] = split[i];
      hisum[i] = sum[i];
    }
  }
  for (i = 1; i < n; i++)
    if (lo[i] < lo[i-1]) {
      lo[i] = lo[i-1];
      losum[i] = losum[i-1];
    }
  for (i = n-1; i > 0; i--)
    if (hi[i] > hi[i+1]) {
      hi[i] = hi[i+1];
      hisum[i] = hisum[i+1];
    }
  int change = 0;
  for (i = 1; i < n; i++)
    if (sum[i] != target[i]) {
      change = 1;
      if (rho == 0) split[i] = 0.5 * (lo[i]+hi[i]);
      else {
        fraction = 1.0*(target[i]-losum[i]) / (hisum[i]-losum[i]);
        split[i] = lo[i] + fraction * (hi[i]-lo[i]);
      }
    }
  return change;
}
double Balance::imbalance_splits()
{
  double *xsplit = comm->xsplit;
  double *ysplit = comm->ysplit;
  double *zsplit = comm->zsplit;
  int nx = comm->procgrid[0];
  int ny = comm->procgrid[1];
  int nz = comm->procgrid[2];
  for (int i = 0; i < nprocs; i++) proccost[i] = 0.0;
  double **x = atom->x;
  int nlocal = atom->nlocal;
  int ix,iy,iz;
  for (int i = 0; i < nlocal; i++) {
    ix = utils::binary_search(x[i][0],nx,xsplit);
    iy = utils::binary_search(x[i][1],ny,ysplit);
    iz = utils::binary_search(x[i][2],nz,zsplit);
    proccost[iz*nx*ny + iy*nx + ix] += 1.0;
  }
  MPI_Allreduce(proccost,allproccost,nprocs,MPI_DOUBLE,MPI_SUM,world);
  double maxcost = 0.0;
  double totalcost = 0.0;
  for (int i = 0; i < nprocs; i++) {
    maxcost = MAX(maxcost,allproccost[i]);
    totalcost += allproccost[i];
  }
  double imbalance = 1.0;
  if (maxcost > 0.0) imbalance = maxcost / (totalcost/nprocs);
  return imbalance;
}
void Balance::dumpout(bigint tstep)
{
  int dimension = domain->dimension;
  int triclinic = domain->triclinic;
  double *lo,*hi;
  if (triclinic == 0) {
    lo = domain->sublo;
    hi = domain->subhi;
  } else {
    lo = domain->sublo_lamda;
    hi = domain->subhi_lamda;
  }
  double box[6];
  box[0] = lo[0]; box[1] = lo[1]; box[2] = lo[2];
  box[3] = hi[0]; box[4] = hi[1]; box[5] = hi[2];
  double **boxall;
  memory->create(boxall,nprocs,6,"balance:dumpout");
  MPI_Allgather(box,6,MPI_DOUBLE,&boxall[0][0],6,MPI_DOUBLE,world);
  if (me) {
    memory->destroy(boxall);
    return;
  }
  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  fmt::print(fp,"ITEM: TIMESTEP\n{}\n",tstep);
  fprintf(fp,"ITEM: NUMBER OF NODES\n");
  if (dimension == 2) fprintf(fp,"%d\n",4*nprocs);
  else fprintf(fp,"%d\n",8*nprocs);
  fprintf(fp,"ITEM: BOX BOUNDS\n");
  fprintf(fp,"%g %g\n",boxlo[0],boxhi[0]);
  fprintf(fp,"%g %g\n",boxlo[1],boxhi[1]);
  fprintf(fp,"%g %g\n",boxlo[2],boxhi[2]);
  fprintf(fp,"ITEM: NODES\n");
  if (triclinic == 0) {
    if (dimension == 2) {
      int m = 0;
      for (int i = 0; i < nprocs; i++) {
        fprintf(fp,"%d %d %g %g %g\n",m+1,1,boxall[i][0],boxall[i][1],0.0);
        fprintf(fp,"%d %d %g %g %g\n",m+2,1,boxall[i][3],boxall[i][1],0.0);
        fprintf(fp,"%d %d %g %g %g\n",m+3,1,boxall[i][3],boxall[i][4],0.0);
        fprintf(fp,"%d %d %g %g %g\n",m+4,1,boxall[i][0],boxall[i][4],0.0);
        m += 4;
      }
    } else {
      int m = 0;
      for (int i = 0; i < nprocs; i++) {
        fprintf(fp,"%d %d %g %g %g\n",m+1,1,
                boxall[i][0],boxall[i][1],boxall[i][2]);
        fprintf(fp,"%d %d %g %g %g\n",m+2,1,
                boxall[i][3],boxall[i][1],boxall[i][2]);
        fprintf(fp,"%d %d %g %g %g\n",m+3,1,
                boxall[i][3],boxall[i][4],boxall[i][2]);
        fprintf(fp,"%d %d %g %g %g\n",m+4,1,
                boxall[i][0],boxall[i][4],boxall[i][2]);
        fprintf(fp,"%d %d %g %g %g\n",m+5,1,
                boxall[i][0],boxall[i][1],boxall[i][5]);
        fprintf(fp,"%d %d %g %g %g\n",m+6,1,
                boxall[i][3],boxall[i][1],boxall[i][5]);
        fprintf(fp,"%d %d %g %g %g\n",m+7,1,
                boxall[i][3],boxall[i][4],boxall[i][5]);
        fprintf(fp,"%d %d %g %g %g\n",m+8,1,
                boxall[i][0],boxall[i][4],boxall[i][5]);
        m += 8;
      }
    }
  } else {
    double (*bc)[3] = domain->corners;
    if (dimension == 2) {
      int m = 0;
      for (int i = 0; i < nprocs; i++) {
        domain->lamda_box_corners(&boxall[i][0],&boxall[i][3]);
        fprintf(fp,"%d %d %g %g %g\n",m+1,1,bc[0][0],bc[0][1],0.0);
        fprintf(fp,"%d %d %g %g %g\n",m+2,1,bc[1][0],bc[1][1],0.0);
        fprintf(fp,"%d %d %g %g %g\n",m+3,1,bc[2][0],bc[2][1],0.0);
        fprintf(fp,"%d %d %g %g %g\n",m+4,1,bc[3][0],bc[3][1],0.0);
        m += 4;
      }
    } else {
      int m = 0;
      for (int i = 0; i < nprocs; i++) {
        domain->lamda_box_corners(&boxall[i][0],&boxall[i][3]);
        fprintf(fp,"%d %d %g %g %g\n",m+1,1,bc[0][0],bc[0][1],bc[0][2]);
        fprintf(fp,"%d %d %g %g %g\n",m+2,1,bc[1][0],bc[1][1],bc[1][2]);
        fprintf(fp,"%d %d %g %g %g\n",m+3,1,bc[2][0],bc[2][1],bc[2][2]);
        fprintf(fp,"%d %d %g %g %g\n",m+4,1,bc[3][0],bc[3][1],bc[3][2]);
        fprintf(fp,"%d %d %g %g %g\n",m+5,1,bc[4][0],bc[4][1],bc[4][2]);
        fprintf(fp,"%d %d %g %g %g\n",m+6,1,bc[5][0],bc[5][1],bc[5][2]);
        fprintf(fp,"%d %d %g %g %g\n",m+7,1,bc[6][0],bc[6][1],bc[6][2]);
        fprintf(fp,"%d %d %g %g %g\n",m+8,1,bc[7][0],bc[7][1],bc[7][2]);
        m += 8;
      }
    }
  }
  fmt::print(fp,"ITEM: TIMESTEP\n{}\n",tstep);
  if (dimension == 2) fprintf(fp,"ITEM: NUMBER OF SQUARES\n");
  else fprintf(fp,"ITEM: NUMBER OF CUBES\n");
  fprintf(fp,"%d\n",nprocs);
  if (dimension == 2) fprintf(fp,"ITEM: SQUARES\n");
  else fprintf(fp,"ITEM: CUBES\n");
  if (dimension == 2) {
    int m = 0;
    for (int i = 0; i < nprocs; i++) {
      fprintf(fp,"%d %d %d %d %d %d\n",i+1,1,m+1,m+2,m+3,m+4);
      m += 4;
    }
  } else {
    int m = 0;
    for (int i = 0; i < nprocs; i++) {
      fprintf(fp,"%d %d %d %d %d %d %d %d %d %d\n",
              i+1,1,m+1,m+2,m+3,m+4,m+5,m+6,m+7,m+8);
      m += 8;
    }
  }
  memory->destroy(boxall);
}
#ifdef BALANCE_DEBUG
void Balance::debug_shift_output(int idim, int m, int np, double *split)
{
  int i;
  const char *dim = nullptr;
  double *boxlo = domain->boxlo;
  double *prd = domain->prd;
  if (bdim[idim] == X) dim = "X";
  else if (bdim[idim] == Y) dim = "Y";
  else if (bdim[idim] == Z) dim = "Z";
  fprintf(stderr,"Dimension %s, Iteration %d\n",dim,m);
  fprintf(stderr,"  Count:");
  for (i = 0; i <= np; i++) fmt::print(stderr," {}",count[i]);
  fprintf(stderr,"\n");
  fprintf(stderr,"  Sum:");
  for (i = 0; i <= np; i++) fmt::print(stderr," {}",sum[i]);
  fprintf(stderr,"\n");
  fprintf(stderr,"  Target:");
  for (i = 0; i <= np; i++) fmt::print(stderr," {}",target[i]);
  fprintf(stderr,"\n");
  fprintf(stderr,"  Actual cut:");
  for (i = 0; i <= np; i++)
    fprintf(stderr," %g",boxlo[bdim[idim]] + split[i]*prd[bdim[idim]]);
  fprintf(stderr,"\n");
  fprintf(stderr,"  Split:");
  for (i = 0; i <= np; i++) fprintf(stderr," %g",split[i]);
  fprintf(stderr,"\n");
  fprintf(stderr,"  Low:");
  for (i = 0; i <= np; i++) fprintf(stderr," %g",lo[i]);
  fprintf(stderr,"\n");
  fprintf(stderr,"  Low-sum:");
  for (i = 0; i <= np; i++) fmt::print(stderr," {}",losum[i]);
  fprintf(stderr,"\n");
  fprintf(stderr,"  Hi:");
  for (i = 0; i <= np; i++) fprintf(stderr," %g",hi[i]);
  fprintf(stderr,"\n");
  fprintf(stderr,"  Hi-sum:");
  for (i = 0; i <= np; i++) fmt::print(stderr," {}",hisum[i]);
  fprintf(stderr,"\n");
  fprintf(stderr,"  Delta:");
  for (i = 0; i < np; i++) fprintf(stderr," %g",split[i+1]-split[i]);
  fprintf(stderr,"\n");
}
#endif
