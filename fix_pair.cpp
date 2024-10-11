#include "fix_pair.h"
#include "atom.h"
#include "error.h"
#include "force.h"
#include "fix.h"
#include "memory.h"
#include "pair.h"
#include "update.h"
#include "fmt/format.h"
using namespace LAMMPS_NS;
using namespace FixConst;
FixPair::FixPair(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) utils::missing_cmd_args(FLERR, "fix pair", error);
  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery < 1) error->all(FLERR,"Illegal fix pair every value: {}", nevery);
  pairname = utils::strdup(arg[4]);
  query_pstyle(lmp);
  if (pstyle == nullptr) error->all(FLERR,"Pair style {} for fix pair not found", pairname);
  nfield = (narg-5) / 2;
  fieldname = new char*[nfield];
  trigger = new int[nfield];
  nfield = 0;
  int iarg = 5;
  while (iarg < narg) {
    if (iarg+2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix pair {}", arg[iarg]), error);
    fieldname[nfield] = utils::strdup(arg[iarg]);
    int flag = utils::inumeric(FLERR,arg[iarg+1],true,lmp);
    if (flag == 0) trigger[nfield] = 0;
    else if (flag == 1) trigger[nfield] = 1;
    else error->all(FLERR,"Illegal fix pair {} command flag: {}", arg[iarg], arg[iarg+1]);
    nfield++;
    iarg += 2;
  }
  triggername = new char*[nfield];
  for (int ifield = 0; ifield < nfield; ifield++) {
    if (trigger[ifield] == 0) triggername[ifield] = nullptr;
    else triggername[ifield] = utils::strdup(fmt::format("{}_flag", fieldname[ifield]));
  }
  triggerptr = new int*[nfield];
  ncols = 0;
  for (int ifield = 0; ifield < nfield; ifield++) {
    int columns = 0;
    pstyle->extract_peratom(fieldname[ifield],columns);
    if (columns) ncols += columns;
    else ncols++;
    if (trigger[ifield]) {
      int dim;
      triggerptr[ifield] = (int *) pstyle->extract(triggername[ifield],dim);
      if (!triggerptr[ifield])
        error->all(FLERR,"Fix pair pair style cannot extract {}", triggername[ifield]);
      if (dim)
        error->all(FLERR,"Fix pair pair style {} trigger {} is not a scalar",
                   pairname, triggername[ifield]);
    }
  }
  peratom_flag = 1;
  if (ncols == 1) size_peratom_cols = 0;
  else size_peratom_cols = ncols;
  peratom_freq = 1;
  vector = nullptr;
  array = nullptr;
  grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);
  int nlocal = atom->nlocal;
  if (ncols == 1) {
    for (int i = 0; i < nlocal; i++)
      vector[i] = 0.0;
  } else {
    for (int i = 0; i < nlocal; i++)
      for (int m = 0; m < ncols; m++)
        array[i][m] = 0.0;
  }
  lasttime = -1;
}
void FixPair::query_pstyle(LAMMPS *lmp) {
    char *cptr=nullptr;
    int nsub = 0;
    if ((cptr = strchr(pairname, ':'))) {
        *cptr = '\0';
        nsub = utils::inumeric(FLERR,cptr+1,false,lmp);
    }
    pstyle = nullptr;
    if (pstyle == nullptr) pstyle = force->pair_match(pairname, 1, nsub);
}
FixPair::~FixPair()
{
  atom->delete_callback(id,Atom::GROW);
  delete[] pairname;
  for (int ifield = 0; ifield < nfield; ifield++) {
    delete[] fieldname[ifield];
    delete[] triggername[ifield];
  }
  delete[] fieldname;
  delete[] trigger;
  delete[] triggername;
  delete[] triggerptr;
  if (ncols == 1) memory->destroy(vector);
  else memory->destroy(array);
}
int FixPair::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= MIN_PRE_FORCE;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  return mask;
}
void FixPair::init()
{
  query_pstyle(lmp);
  if (pstyle == nullptr) error->all(FLERR,"Pair style {} for fix pair not found", pairname);
}
void FixPair::setup(int vflag)
{
  post_force(vflag);
}
void FixPair::min_setup(int vflag)
{
  setup(vflag);
}
void FixPair::setup_pre_force(int vflag)
{
  pre_force(vflag);
}
void FixPair::pre_force(int )
{
  if (update->ntimestep % nevery) return;
  if (update->ntimestep == lasttime) return;
  for (int ifield = 0; ifield < nfield; ifield++)
    if (trigger[ifield]) *(triggerptr[ifield]) = 1;
}
void FixPair::min_pre_force(int vflag)
{
  pre_force(vflag);
}
void FixPair::post_force(int )
{
  if (update->ntimestep % nevery) return;
  if (update->ntimestep == lasttime) return;
  lasttime = update->ntimestep;
  const int nlocal = atom->nlocal;
  int icol = 0;
  int columns;
  for (int ifield = 0; ifield < nfield; ifield++) {
    void *pvoid = pstyle->extract_peratom(fieldname[ifield],columns);
    if ((pvoid == nullptr) && (nlocal > 0))
      error->one(FLERR, "Fix pair cannot extract property {} from pair style", fieldname[ifield]);
    if (columns == 0) {
      double *pvector = (double *) pvoid;
      if (ncols == 1) {
        for (int i = 0; i < nlocal; i++)
          vector[i] = pvector[i];
      } else {
        for (int i = 0; i < nlocal; i++)
          array[i][icol] = pvector[i];
      }
      icol++;
    } else {
      double **parray = (double **) pvoid;
      for (int i = 0; i < nlocal; i++)
        for (int m = 0; m < columns; m++) {
          array[i][icol] = parray[i][m];
          icol++;
        }
    }
  }
  for (int ifield = 0; ifield < nfield; ifield++)
    if (trigger[ifield]) *(triggerptr[ifield]) = 0;
}
void FixPair::min_post_force(int vflag)
{
  post_force(vflag);
}
void FixPair::grow_arrays(int nmax)
{
  if (ncols == 1) {
    memory->grow(vector,nmax,"store/state:vector");
    vector_atom = vector;
  } else {
    memory->grow(array,nmax,ncols,"store/state:array");
    array_atom = array;
  }
}
void FixPair::copy_arrays(int i, int j, int )
{
  if (ncols == 1) {
    vector[j] = vector[i];
  } else {
    for (int m = 0; m < ncols; m++)
      array[j][m] = array[i][m];
  }
}
int FixPair::pack_exchange(int i, double *buf)
{
  if (ncols == 1) {
    buf[0] = vector[i];
  } else {
    for (int m = 0; m < ncols; m++)
      buf[m] = array[i][m];
  }
  return ncols;
}
int FixPair::unpack_exchange(int nlocal, double *buf)
{
  if (ncols == 1) {
    vector[nlocal] = buf[0];
  } else {
    for (int m = 0; m < ncols; m++)
      array[nlocal][m] = buf[m];
  }
  return ncols;
}
double FixPair::memory_usage()
{
  double bytes = 0.0;
  if (ncols == 1) bytes += (double)atom->nmax * sizeof(double);
  else bytes += (double)atom->nmax*ncols * sizeof(double);
  return bytes;
}
