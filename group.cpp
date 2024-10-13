#include "group.h"
#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "input.h"
#include "math_extra.h"
#include "math_eigen.h"
#include "memory.h"
#include "modify.h"
#include "region.h"
#include "tokenizer.h"
#include "exceptions.h"
#include <cmath>
#include <cstring>
#include <map>
#include <utility>
using namespace LAMMPS_NS;
static constexpr int MAX_GROUP = 32;
static constexpr double EPSILON = 1.0e-6;
enum{NONE,TYPE,MOLECULE,ID};
enum{LT,LE,GT,GE,EQ,NEQ,BETWEEN};
#define BIG 1.0e20
Group::Group(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  names = new char*[MAX_GROUP];
  bitmask = new int[MAX_GROUP];
  inversemask = new int[MAX_GROUP];
  dynamic = new int[MAX_GROUP];
  for (int i = 0; i < MAX_GROUP; i++) names[i] = nullptr;
  for (int i = 0; i < MAX_GROUP; i++) bitmask[i] = 1 << i;
  for (int i = 0; i < MAX_GROUP; i++) inversemask[i] = bitmask[i] ^ ~0;
  for (int i = 0; i < MAX_GROUP; i++) dynamic[i] = 0;
  names[0] = utils::strdup("all");
  ngroup = 1;
}
Group::~Group()
{
  for (int i = 0; i < MAX_GROUP; i++) delete [] names[i];
  delete [] names;
  delete [] bitmask;
  delete [] inversemask;
  delete [] dynamic;
}
int Group::find(const std::string &name)
{
  for (int igroup = 0; igroup < MAX_GROUP; igroup++)
    if (names[igroup] && (name == names[igroup])) return igroup;
  return -1;
}
int Group::find_unused()
{
  for (int igroup = 0; igroup < MAX_GROUP; igroup++)
    if (names[igroup] == nullptr) return igroup;
  return -1;
}
bigint Group::count_all()
{
  bigint nme = atom->nlocal;
  bigint nall;
  MPI_Allreduce(&nme,&nall,1,MPI_LMP_BIGINT,MPI_SUM,world);
  return nall;
}
bigint Group::count(int igroup)
{
  int groupbit = bitmask[igroup];
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int n = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) n++;
  bigint nsingle = n;
  bigint nall;
  MPI_Allreduce(&nsingle,&nall,1,MPI_LMP_BIGINT,MPI_SUM,world);
  return nall;
}
bigint Group::count(int igroup, Region *region)
{
  region->prematch();
  const int groupbit = bitmask[igroup];
  double **x = atom->x;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;
  int n = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) n++;
  bigint nsingle = n;
  bigint nall;
  MPI_Allreduce(&nsingle,&nall,1,MPI_LMP_BIGINT,MPI_SUM,world);
  return nall;
}
double Group::mass(int igroup)
{
  int groupbit = bitmask[igroup];
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double one = 0.0;
  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) one += rmass[i];
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) one += mass[type[i]];
  }
  double all;
  MPI_Allreduce(&one,&all,1,MPI_DOUBLE,MPI_SUM,world);
  return all;
}
double Group::mass(int igroup, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();
  double **x = atom->x;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double one = 0.0;
  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
        one += rmass[i];
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
        one += mass[type[i]];
  }
  double all;
  MPI_Allreduce(&one,&all,1,MPI_DOUBLE,MPI_SUM,world);
  return all;
}
double Group::charge(int igroup)
{
  int groupbit = bitmask[igroup];
  double *q = atom->q;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double qone = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) qone += q[i];
  double qall;
  MPI_Allreduce(&qone,&qall,1,MPI_DOUBLE,MPI_SUM,world);
  return qall;
}
double Group::charge(int igroup, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();
  double **x = atom->x;
  double *q = atom->q;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double qone = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
      qone += q[i];
  double qall;
  MPI_Allreduce(&qone,&qall,1,MPI_DOUBLE,MPI_SUM,world);
  return qall;
}
void Group::bounds(int igroup, double *minmax)
{
  int groupbit = bitmask[igroup];
  double extent[6];
  extent[0] = extent[2] = extent[4] = BIG;
  extent[1] = extent[3] = extent[5] = -BIG;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      extent[0] = MIN(extent[0],x[i][0]);
      extent[1] = MAX(extent[1],x[i][0]);
      extent[2] = MIN(extent[2],x[i][1]);
      extent[3] = MAX(extent[3],x[i][1]);
      extent[4] = MIN(extent[4],x[i][2]);
      extent[5] = MAX(extent[5],x[i][2]);
    }
  }
  extent[0] = -extent[0];
  extent[2] = -extent[2];
  extent[4] = -extent[4];
  MPI_Allreduce(extent,minmax,6,MPI_DOUBLE,MPI_MAX,world);
  minmax[0] = -minmax[0];
  minmax[2] = -minmax[2];
  minmax[4] = -minmax[4];
}
void Group::bounds(int igroup, double *minmax, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();
  double extent[6];
  extent[0] = extent[2] = extent[4] = BIG;
  extent[1] = extent[3] = extent[5] = -BIG;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
      extent[0] = MIN(extent[0],x[i][0]);
      extent[1] = MAX(extent[1],x[i][0]);
      extent[2] = MIN(extent[2],x[i][1]);
      extent[3] = MAX(extent[3],x[i][1]);
      extent[4] = MIN(extent[4],x[i][2]);
      extent[5] = MAX(extent[5],x[i][2]);
    }
  }
  extent[0] = -extent[0];
  extent[2] = -extent[2];
  extent[4] = -extent[4];
  MPI_Allreduce(extent,minmax,6,MPI_DOUBLE,MPI_MAX,world);
  minmax[0] = -minmax[0];
  minmax[2] = -minmax[2];
  minmax[4] = -minmax[4];
}
void Group::xcm(int igroup, double masstotal, double *cm)
{
  int groupbit = bitmask[igroup];
  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double cmone[3];
  cmone[0] = cmone[1] = cmone[2] = 0.0;
  double massone;
  double unwrap[3];
  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        massone = rmass[i];
        domain->unmap(x[i],image[i],unwrap);
        cmone[0] += unwrap[0] * massone;
        cmone[1] += unwrap[1] * massone;
        cmone[2] += unwrap[2] * massone;
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        massone = mass[type[i]];
        domain->unmap(x[i],image[i],unwrap);
        cmone[0] += unwrap[0] * massone;
        cmone[1] += unwrap[1] * massone;
        cmone[2] += unwrap[2] * massone;
      }
  }
  MPI_Allreduce(cmone,cm,3,MPI_DOUBLE,MPI_SUM,world);
  if (masstotal > 0.0) {
    cm[0] /= masstotal;
    cm[1] /= masstotal;
    cm[2] /= masstotal;
  }
}
void Group::xcm(int igroup, double masstotal, double *cm, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();
  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double cmone[3];
  cmone[0] = cmone[1] = cmone[2] = 0.0;
  double massone;
  double unwrap[3];
  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
        massone = rmass[i];
        domain->unmap(x[i],image[i],unwrap);
        cmone[0] += unwrap[0] * massone;
        cmone[1] += unwrap[1] * massone;
        cmone[2] += unwrap[2] * massone;
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
        massone = mass[type[i]];
        domain->unmap(x[i],image[i],unwrap);
        cmone[0] += unwrap[0] * massone;
        cmone[1] += unwrap[1] * massone;
        cmone[2] += unwrap[2] * massone;
      }
  }
  MPI_Allreduce(cmone,cm,3,MPI_DOUBLE,MPI_SUM,world);
  if (masstotal > 0.0) {
    cm[0] /= masstotal;
    cm[1] /= masstotal;
    cm[2] /= masstotal;
  }
}
void Group::vcm(int igroup, double masstotal, double *cm)
{
  int groupbit = bitmask[igroup];
  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double p[3],massone;
  p[0] = p[1] = p[2] = 0.0;
  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        massone = rmass[i];
        p[0] += v[i][0]*massone;
        p[1] += v[i][1]*massone;
        p[2] += v[i][2]*massone;
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        massone = mass[type[i]];
        p[0] += v[i][0]*massone;
        p[1] += v[i][1]*massone;
        p[2] += v[i][2]*massone;
      }
  }
  MPI_Allreduce(p,cm,3,MPI_DOUBLE,MPI_SUM,world);
  if (masstotal > 0.0) {
    cm[0] /= masstotal;
    cm[1] /= masstotal;
    cm[2] /= masstotal;
  }
}
void Group::vcm(int igroup, double masstotal, double *cm, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double p[3],massone;
  p[0] = p[1] = p[2] = 0.0;
  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
        massone = rmass[i];
        p[0] += v[i][0]*massone;
        p[1] += v[i][1]*massone;
        p[2] += v[i][2]*massone;
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
        massone = mass[type[i]];
        p[0] += v[i][0]*massone;
        p[1] += v[i][1]*massone;
        p[2] += v[i][2]*massone;
      }
  }
  MPI_Allreduce(p,cm,3,MPI_DOUBLE,MPI_SUM,world);
  if (masstotal > 0.0) {
    cm[0] /= masstotal;
    cm[1] /= masstotal;
    cm[2] /= masstotal;
  }
}
void Group::fcm(int igroup, double *cm)
{
  int groupbit = bitmask[igroup];
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double flocal[3];
  flocal[0] = flocal[1] = flocal[2] = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      flocal[0] += f[i][0];
      flocal[1] += f[i][1];
      flocal[2] += f[i][2];
    }
  MPI_Allreduce(flocal,cm,3,MPI_DOUBLE,MPI_SUM,world);
}
void Group::fcm(int igroup, double *cm, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double flocal[3];
  flocal[0] = flocal[1] = flocal[2] = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
      flocal[0] += f[i][0];
      flocal[1] += f[i][1];
      flocal[2] += f[i][2];
    }
  MPI_Allreduce(flocal,cm,3,MPI_DOUBLE,MPI_SUM,world);
}
double Group::ke(int igroup)
{
  int groupbit = bitmask[igroup];
  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double one = 0.0;
  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        one += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) *
          rmass[i];
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        one += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) *
          mass[type[i]];
  }
  double all;
  MPI_Allreduce(&one,&all,1,MPI_DOUBLE,MPI_SUM,world);
  all *= 0.5 * force->mvv2e;
  return all;
}
double Group::ke(int igroup, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double one = 0.0;
  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
        one += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) *
          rmass[i];
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
        one += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) *
          mass[type[i]];
  }
  double all;
  MPI_Allreduce(&one,&all,1,MPI_DOUBLE,MPI_SUM,world);
  all *= 0.5 * force->mvv2e;
  return all;
}
double Group::gyration(int igroup, double masstotal, double *cm)
{
  int groupbit = bitmask[igroup];
  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double dx,dy,dz,massone;
  double unwrap[3];
  double rg = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - cm[0];
      dy = unwrap[1] - cm[1];
      dz = unwrap[2] - cm[2];
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      rg += (dx*dx + dy*dy + dz*dz) * massone;
    }
  double rg_all;
  MPI_Allreduce(&rg,&rg_all,1,MPI_DOUBLE,MPI_SUM,world);
  if (masstotal > 0.0) return sqrt(rg_all/masstotal);
  return 0.0;
}
double Group::gyration(int igroup, double masstotal, double *cm, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();
  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double dx,dy,dz,massone;
  double unwrap[3];
  double rg = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - cm[0];
      dy = unwrap[1] - cm[1];
      dz = unwrap[2] - cm[2];
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      rg += (dx*dx + dy*dy + dz*dz) * massone;
    }
  double rg_all;
  MPI_Allreduce(&rg,&rg_all,1,MPI_DOUBLE,MPI_SUM,world);
  if (masstotal > 0.0) return sqrt(rg_all/masstotal);
  return 0.0;
}
void Group::angmom(int igroup, double *cm, double *lmom)
{
  int groupbit = bitmask[igroup];
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double dx,dy,dz,massone;
  double unwrap[3];
  double p[3];
  p[0] = p[1] = p[2] = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - cm[0];
      dy = unwrap[1] - cm[1];
      dz = unwrap[2] - cm[2];
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      p[0] += massone * (dy*v[i][2] - dz*v[i][1]);
      p[1] += massone * (dz*v[i][0] - dx*v[i][2]);
      p[2] += massone * (dx*v[i][1] - dy*v[i][0]);
    }
  MPI_Allreduce(p,lmom,3,MPI_DOUBLE,MPI_SUM,world);
}
void Group::angmom(int igroup, double *cm, double *lmom, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double dx,dy,dz,massone;
  double unwrap[3];
  double p[3];
  p[0] = p[1] = p[2] = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - cm[0];
      dy = unwrap[1] - cm[1];
      dz = unwrap[2] - cm[2];
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      p[0] += massone * (dy*v[i][2] - dz*v[i][1]);
      p[1] += massone * (dz*v[i][0] - dx*v[i][2]);
      p[2] += massone * (dx*v[i][1] - dy*v[i][0]);
    }
  MPI_Allreduce(p,lmom,3,MPI_DOUBLE,MPI_SUM,world);
}
void Group::torque(int igroup, double *cm, double *tq)
{
  int groupbit = bitmask[igroup];
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;
  double dx,dy,dz;
  double unwrap[3];
  double tlocal[3];
  tlocal[0] = tlocal[1] = tlocal[2] = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - cm[0];
      dy = unwrap[1] - cm[1];
      dz = unwrap[2] - cm[2];
      tlocal[0] += dy*f[i][2] - dz*f[i][1];
      tlocal[1] += dz*f[i][0] - dx*f[i][2];
      tlocal[2] += dx*f[i][1] - dy*f[i][0];
    }
  MPI_Allreduce(tlocal,tq,3,MPI_DOUBLE,MPI_SUM,world);
}
void Group::torque(int igroup, double *cm, double *tq, Region *region)
{
  int groupbit = bitmask[igroup];
  region->prematch();
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;
  double dx,dy,dz;
  double unwrap[3];
  double tlocal[3];
  tlocal[0] = tlocal[1] = tlocal[2] = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - cm[0];
      dy = unwrap[1] - cm[1];
      dz = unwrap[2] - cm[2];
      tlocal[0] += dy*f[i][2] - dz*f[i][1];
      tlocal[1] += dz*f[i][0] - dx*f[i][2];
      tlocal[2] += dx*f[i][1] - dy*f[i][0];
    }
  MPI_Allreduce(tlocal,tq,3,MPI_DOUBLE,MPI_SUM,world);
}
void Group::inertia(int igroup, double *cm, double itensor[3][3])
{
  int i,j;
  int groupbit = bitmask[igroup];
  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double dx,dy,dz,massone;
  double unwrap[3];
  double ione[3][3];
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      ione[i][j] = 0.0;
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - cm[0];
      dy = unwrap[1] - cm[1];
      dz = unwrap[2] - cm[2];
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      ione[0][0] += massone * (dy*dy + dz*dz);
      ione[1][1] += massone * (dx*dx + dz*dz);
      ione[2][2] += massone * (dx*dx + dy*dy);
      ione[0][1] -= massone * dx*dy;
      ione[1][2] -= massone * dy*dz;
      ione[0][2] -= massone * dx*dz;
    }
  ione[1][0] = ione[0][1];
  ione[2][1] = ione[1][2];
  ione[2][0] = ione[0][2];
  MPI_Allreduce(&ione[0][0],&itensor[0][0],9,MPI_DOUBLE,MPI_SUM,world);
}
void Group::inertia(int igroup, double *cm, double itensor[3][3], Region *region)
{
  int i,j;
  int groupbit = bitmask[igroup];
  region->prematch();
  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double dx,dy,dz,massone;
  double unwrap[3];
  double ione[3][3];
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      ione[i][j] = 0.0;
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - cm[0];
      dy = unwrap[1] - cm[1];
      dz = unwrap[2] - cm[2];
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      ione[0][0] += massone * (dy*dy + dz*dz);
      ione[1][1] += massone * (dx*dx + dz*dz);
      ione[2][2] += massone * (dx*dx + dy*dy);
      ione[0][1] -= massone * dx*dy;
      ione[1][2] -= massone * dy*dz;
      ione[0][2] -= massone * dx*dz;
    }
  ione[1][0] = ione[0][1];
  ione[2][1] = ione[1][2];
  ione[2][0] = ione[0][2];
  MPI_Allreduce(&ione[0][0],&itensor[0][0],9,MPI_DOUBLE,MPI_SUM,world);
}
void Group::omega(double *angmom, double inertia[3][3], double *w)
{
  double idiag[3],ex[3],ey[3],ez[3],cross[3];
  double evectors[3][3],inverse[3][3];
  double determinant = inertia[0][0]*inertia[1][1]*inertia[2][2] +
    inertia[0][1]*inertia[1][2]*inertia[2][0] +
    inertia[0][2]*inertia[1][0]*inertia[2][1] -
    inertia[0][0]*inertia[1][2]*inertia[2][1] -
    inertia[0][1]*inertia[1][0]*inertia[2][2] -
    inertia[2][0]*inertia[1][1]*inertia[0][2];
  if (determinant > EPSILON) {
    inverse[0][0] = inertia[1][1]*inertia[2][2] - inertia[1][2]*inertia[2][1];
    inverse[0][1] = -(inertia[0][1]*inertia[2][2] -
                      inertia[0][2]*inertia[2][1]);
    inverse[0][2] = inertia[0][1]*inertia[1][2] - inertia[0][2]*inertia[1][1];
    inverse[1][0] = -(inertia[1][0]*inertia[2][2] -
                      inertia[1][2]*inertia[2][0]);
    inverse[1][1] = inertia[0][0]*inertia[2][2] - inertia[0][2]*inertia[2][0];
    inverse[1][2] = -(inertia[0][0]*inertia[1][2] -
                      inertia[0][2]*inertia[1][0]);
    inverse[2][0] = inertia[1][0]*inertia[2][1] - inertia[1][1]*inertia[2][0];
    inverse[2][1] = -(inertia[0][0]*inertia[2][1] -
                      inertia[0][1]*inertia[2][0]);
    inverse[2][2] = inertia[0][0]*inertia[1][1] - inertia[0][1]*inertia[1][0];
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        inverse[i][j] /= determinant;
    w[0] = inverse[0][0]*angmom[0] + inverse[0][1]*angmom[1] +
      inverse[0][2]*angmom[2];
    w[1] = inverse[1][0]*angmom[0] + inverse[1][1]*angmom[1] +
      inverse[1][2]*angmom[2];
    w[2] = inverse[2][0]*angmom[0] + inverse[2][1]*angmom[1] +
      inverse[2][2]*angmom[2];
  } else {
    int ierror = MathEigen::jacobi3(inertia, idiag, evectors);
    if (ierror) error->all(FLERR, "Insufficient Jacobi rotations for group::omega");
    ex[0] = evectors[0][0];
    ex[1] = evectors[1][0];
    ex[2] = evectors[2][0];
    ey[0] = evectors[0][1];
    ey[1] = evectors[1][1];
    ey[2] = evectors[2][1];
    ez[0] = evectors[0][2];
    ez[1] = evectors[1][2];
    ez[2] = evectors[2][2];
    MathExtra::cross3(ex,ey,cross);
    if (MathExtra::dot3(cross,ez) < 0.0) MathExtra::negate3(ez);
    double max;
    max = MAX(idiag[0],idiag[1]);
    max = MAX(max,idiag[2]);
    if (idiag[0] < EPSILON*max) idiag[0] = 0.0;
    if (idiag[1] < EPSILON*max) idiag[1] = 0.0;
    if (idiag[2] < EPSILON*max) idiag[2] = 0.0;
    MathExtra::angmom_to_omega(angmom,ex,ey,ez,idiag,w);
  }
}
