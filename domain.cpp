#include "domain.h"
#include "style_region.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "lattice.h"
#include "memory.h"
#include "modify.h"
#include "region.h"
#include "universe.h"
#include "update.h"
#include <cstring>
#include <cmath>
using namespace LAMMPS_NS;
#define BIG 1.0e20
#define SMALL 1.0e-4
#define DELTAREGION 4
#define BONDSTRETCH 1.1
template <typename T> static Region *region_creator(LAMMPS *lmp, int narg, char ** arg)
{
  return new T(lmp, narg, arg);
}
Domain::Domain(LAMMPS *lmp) : Pointers(lmp)
{
  box_exist = 0;
  box_change = 0;
  deform_flag = deform_vremap = deform_groupbit = 0;
  dimension = 3;
  nonperiodic = 0;
  xperiodic = yperiodic = zperiodic = 1;
  periodicity[0] = xperiodic;
  periodicity[1] = yperiodic;
  periodicity[2] = zperiodic;
  boundary[0][0] = boundary[0][1] = 0;
  boundary[1][0] = boundary[1][1] = 0;
  boundary[2][0] = boundary[2][1] = 0;
  minxlo = minxhi = 0.0;
  minylo = minyhi = 0.0;
  minzlo = minzhi = 0.0;
  triclinic = 0;
  boxlo[0] = boxlo[1] = boxlo[2] = -0.5;
  boxhi[0] = boxhi[1] = boxhi[2] = 0.5;
  xy = xz = yz = 0.0;
  h[3] = h[4] = h[5] = 0.0;
  h_inv[3] = h_inv[4] = h_inv[5] = 0.0;
  h_rate[0] = h_rate[1] = h_rate[2] =
    h_rate[3] = h_rate[4] = h_rate[5] = 0.0;
  h_ratelo[0] = h_ratelo[1] = h_ratelo[2] = 0.0;
  prd_lamda[0] = prd_lamda[1] = prd_lamda[2] = 1.0;
  prd_half_lamda[0] = prd_half_lamda[1] = prd_half_lamda[2] = 0.5;
  boxlo_lamda[0] = boxlo_lamda[1] = boxlo_lamda[2] = 0.0;
  boxhi_lamda[0] = boxhi_lamda[1] = boxhi_lamda[2] = 1.0;
  lattice = nullptr;
  auto args = new char*[2];
  args[0] = (char *) "none";
  args[1] = (char *) "1.0";
  set_lattice(2,args);
  delete[] args;
  copymode = 0;
  region_map = new RegionCreatorMap();
#define REGION_CLASS 
#define RegionStyle(key,Class) \
  (*region_map)[#key] = &region_creator<Class>;
#include "style_region.h"
#undef RegionStyle
#undef REGION_CLASS
}
Domain::~Domain()
{
  if (copymode) return;
  for (auto &reg : regions) delete reg;
  regions.clear();
  delete lattice;
  delete region_map;
}
void Domain::init()
{
  box_change_size = box_change_shape = box_change_domain = 0;
  int box_change_x=0, box_change_y=0, box_change_z=0;
  int box_change_yz=0, box_change_xz=0, box_change_xy=0;
  const auto &fixes = modify->get_fix_list();
  if (nonperiodic == 2) box_change_size = 1;
  for (const auto &fix : fixes) {
    if (fix->box_change & Fix::BOX_CHANGE_SIZE) box_change_size = 1;
    if (fix->box_change & Fix::BOX_CHANGE_SHAPE) box_change_shape = 1;
    if (fix->box_change & Fix::BOX_CHANGE_DOMAIN) box_change_domain = 1;
    if (fix->box_change & Fix::BOX_CHANGE_X) box_change_x++;
    if (fix->box_change & Fix::BOX_CHANGE_Y) box_change_y++;
    if (fix->box_change & Fix::BOX_CHANGE_Z) box_change_z++;
    if (fix->box_change & Fix::BOX_CHANGE_YZ) box_change_yz++;
    if (fix->box_change & Fix::BOX_CHANGE_XZ) box_change_xz++;
    if (fix->box_change & Fix::BOX_CHANGE_XY) box_change_xy++;
  }
  std::string mesg = "Must not have multiple fixes change box parameter ";
#define CHECK_BOX_FIX_ERROR(par) \
  if (box_change_ ## par > 1) error->all(FLERR,(mesg + #par))
  CHECK_BOX_FIX_ERROR(x);
  CHECK_BOX_FIX_ERROR(y);
  CHECK_BOX_FIX_ERROR(z);
  CHECK_BOX_FIX_ERROR(yz);
  CHECK_BOX_FIX_ERROR(xz);
  CHECK_BOX_FIX_ERROR(xy);
#undef CHECK_BOX_FIX_ERROR
  box_change = 0;
  if (box_change_size || box_change_shape || box_change_domain) box_change = 1;
  deform_flag = deform_vremap = deform_groupbit = 0;
  for (auto &reg : regions) reg->init();
}
void Domain::set_initial_box(int expandflag)
{
  if (boxlo[0] >= boxhi[0] || boxlo[1] >= boxhi[1] || boxlo[2] >= boxhi[2])
    error->one(FLERR,"Box bounds are invalid or missing");
  if (dimension == 2 && (xz != 0.0 || yz != 0.0))
    error->all(FLERR,"Cannot skew triclinic box in z for 2d simulation");
  if (triclinic) {
    if ((fabs(xy/(boxhi[1]-boxlo[1])) > 0.5 && yperiodic) ||
        ((fabs(xz)+fabs(yz))/(boxhi[2]-boxlo[2]) > 0.5 && zperiodic)) {
      if (comm->me == 0)
        error->warning(FLERR,"Triclinic box skew is large. LAMMPS will run inefficiently.");
    }
  }
  small[0] = SMALL * (boxhi[0] - boxlo[0]);
  small[1] = SMALL * (boxhi[1] - boxlo[1]);
  small[2] = SMALL * (boxhi[2] - boxlo[2]);
  if (!expandflag) return;
  if (boundary[0][0] == 2) boxlo[0] -= small[0];
  else if (boundary[0][0] == 3) minxlo = boxlo[0];
  if (boundary[0][1] == 2) boxhi[0] += small[0];
  else if (boundary[0][1] == 3) minxhi = boxhi[0];
  if (boundary[1][0] == 2) boxlo[1] -= small[1];
  else if (boundary[1][0] == 3) minylo = boxlo[1];
  if (boundary[1][1] == 2) boxhi[1] += small[1];
  else if (boundary[1][1] == 3) minyhi = boxhi[1];
  if (boundary[2][0] == 2) boxlo[2] -= small[2];
  else if (boundary[2][0] == 3) minzlo = boxlo[2];
  if (boundary[2][1] == 2) boxhi[2] += small[2];
  else if (boundary[2][1] == 3) minzhi = boxhi[2];
}
void Domain::set_global_box()
{
  prd[0] = xprd = boxhi[0] - boxlo[0];
  prd[1] = yprd = boxhi[1] - boxlo[1];
  prd[2] = zprd = boxhi[2] - boxlo[2];
  h[0] = xprd;
  h[1] = yprd;
  h[2] = zprd;
  h_inv[0] = 1.0/h[0];
  h_inv[1] = 1.0/h[1];
  h_inv[2] = 1.0/h[2];
  prd_half[0] = xprd_half = 0.5*xprd;
  prd_half[1] = yprd_half = 0.5*yprd;
  prd_half[2] = zprd_half = 0.5*zprd;
  if (triclinic) {
    h[3] = yz;
    h[4] = xz;
    h[5] = xy;
    h_inv[3] = -h[3] / (h[1]*h[2]);
    h_inv[4] = (h[3]*h[5] - h[1]*h[4]) / (h[0]*h[1]*h[2]);
    h_inv[5] = -h[5] / (h[0]*h[1]);
    boxlo_bound[0] = MIN(boxlo[0],boxlo[0]+xy);
    boxlo_bound[0] = MIN(boxlo_bound[0],boxlo_bound[0]+xz);
    boxlo_bound[1] = MIN(boxlo[1],boxlo[1]+yz);
    boxlo_bound[2] = boxlo[2];
    boxhi_bound[0] = MAX(boxhi[0],boxhi[0]+xy);
    boxhi_bound[0] = MAX(boxhi_bound[0],boxhi_bound[0]+xz);
    boxhi_bound[1] = MAX(boxhi[1],boxhi[1]+yz);
    boxhi_bound[2] = boxhi[2];
  }
}
void Domain::set_lamda_box()
{
  if (comm->layout != Comm::LAYOUT_TILED) {
    int *myloc = comm->myloc;
    double *xsplit = comm->xsplit;
    double *ysplit = comm->ysplit;
    double *zsplit = comm->zsplit;
    sublo_lamda[0] = xsplit[myloc[0]];
    subhi_lamda[0] = xsplit[myloc[0]+1];
    sublo_lamda[1] = ysplit[myloc[1]];
    subhi_lamda[1] = ysplit[myloc[1]+1];
    sublo_lamda[2] = zsplit[myloc[2]];
    subhi_lamda[2] = zsplit[myloc[2]+1];
  } else {
    double (*mysplit)[2] = comm->mysplit;
    sublo_lamda[0] = mysplit[0][0];
    subhi_lamda[0] = mysplit[0][1];
    sublo_lamda[1] = mysplit[1][0];
    subhi_lamda[1] = mysplit[1][1];
    sublo_lamda[2] = mysplit[2][0];
    subhi_lamda[2] = mysplit[2][1];
  }
}
void Domain::set_local_box()
{
  if (triclinic) return;
  if (comm->layout != Comm::LAYOUT_TILED) {
    int *myloc = comm->myloc;
    int *procgrid = comm->procgrid;
    double *xsplit = comm->xsplit;
    double *ysplit = comm->ysplit;
    double *zsplit = comm->zsplit;
    sublo[0] = boxlo[0] + xprd*xsplit[myloc[0]];
    if (myloc[0] < procgrid[0]-1) subhi[0] = boxlo[0] + xprd*xsplit[myloc[0]+1];
    else subhi[0] = boxhi[0];
    sublo[1] = boxlo[1] + yprd*ysplit[myloc[1]];
    if (myloc[1] < procgrid[1]-1) subhi[1] = boxlo[1] + yprd*ysplit[myloc[1]+1];
    else subhi[1] = boxhi[1];
    sublo[2] = boxlo[2] + zprd*zsplit[myloc[2]];
    if (myloc[2] < procgrid[2]-1) subhi[2] = boxlo[2] + zprd*zsplit[myloc[2]+1];
    else subhi[2] = boxhi[2];
  } else {
    double (*mysplit)[2] = comm->mysplit;
    sublo[0] = boxlo[0] + xprd*mysplit[0][0];
    if (mysplit[0][1] < 1.0) subhi[0] = boxlo[0] + xprd*mysplit[0][1];
    else subhi[0] = boxhi[0];
    sublo[1] = boxlo[1] + yprd*mysplit[1][0];
    if (mysplit[1][1] < 1.0) subhi[1] = boxlo[1] + yprd*mysplit[1][1];
    else subhi[1] = boxhi[1];
    sublo[2] = boxlo[2] + zprd*mysplit[2][0];
    if (mysplit[2][1] < 1.0) subhi[2] = boxlo[2] + zprd*mysplit[2][1];
    else subhi[2] = boxhi[2];
  }
}
void Domain::reset_box()
{
  if (atom->natoms == 0) return;
  if (nonperiodic == 2) {
    double extent[3][2],all[3][2];
    extent[2][0] = extent[1][0] = extent[0][0] = BIG;
    extent[2][1] = extent[1][1] = extent[0][1] = -BIG;
    double **x = atom->x;
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++) {
      extent[0][0] = MIN(extent[0][0],x[i][0]);
      extent[0][1] = MAX(extent[0][1],x[i][0]);
      extent[1][0] = MIN(extent[1][0],x[i][1]);
      extent[1][1] = MAX(extent[1][1],x[i][1]);
      extent[2][0] = MIN(extent[2][0],x[i][2]);
      extent[2][1] = MAX(extent[2][1],x[i][2]);
    }
    extent[0][0] = -extent[0][0];
    extent[1][0] = -extent[1][0];
    extent[2][0] = -extent[2][0];
    MPI_Allreduce(extent,all,6,MPI_DOUBLE,MPI_MAX,world);
    if (triclinic) lamda2x(atom->nlocal);
    if (triclinic == 0) {
      if (xperiodic == 0) {
        if (boundary[0][0] == 2) boxlo[0] = -all[0][0] - small[0];
        else if (boundary[0][0] == 3)
          boxlo[0] = MIN(-all[0][0]-small[0],minxlo);
        if (boundary[0][1] == 2) boxhi[0] = all[0][1] + small[0];
        else if (boundary[0][1] == 3) boxhi[0] = MAX(all[0][1]+small[0],minxhi);
        if (boxlo[0] > boxhi[0]) error->all(FLERR,"Illegal simulation box");
      }
      if (yperiodic == 0) {
        if (boundary[1][0] == 2) boxlo[1] = -all[1][0] - small[1];
        else if (boundary[1][0] == 3)
          boxlo[1] = MIN(-all[1][0]-small[1],minylo);
        if (boundary[1][1] == 2) boxhi[1] = all[1][1] + small[1];
        else if (boundary[1][1] == 3) boxhi[1] = MAX(all[1][1]+small[1],minyhi);
        if (boxlo[1] > boxhi[1]) error->all(FLERR,"Illegal simulation box");
      }
      if (zperiodic == 0) {
        if (boundary[2][0] == 2) boxlo[2] = -all[2][0] - small[2];
        else if (boundary[2][0] == 3)
          boxlo[2] = MIN(-all[2][0]-small[2],minzlo);
        if (boundary[2][1] == 2) boxhi[2] = all[2][1] + small[2];
        else if (boundary[2][1] == 3) boxhi[2] = MAX(all[2][1]+small[2],minzhi);
        if (boxlo[2] > boxhi[2]) error->all(FLERR,"Illegal simulation box");
      }
    } else {
      double lo[3],hi[3];
      if (xperiodic == 0) {
        lo[0] = -all[0][0]; lo[1] = 0.0; lo[2] = 0.0;
        lamda2x(lo,lo);
        hi[0] = all[0][1]; hi[1] = 0.0; hi[2] = 0.0;
        lamda2x(hi,hi);
        if (boundary[0][0] == 2) boxlo[0] = lo[0] - small[0];
        else if (boundary[0][0] == 3) boxlo[0] = MIN(lo[0]-small[0],minxlo);
        if (boundary[0][1] == 2) boxhi[0] = hi[0] + small[0];
        else if (boundary[0][1] == 3) boxhi[0] = MAX(hi[0]+small[0],minxhi);
        if (boxlo[0] > boxhi[0]) error->all(FLERR,"Illegal simulation box");
      }
      if (yperiodic == 0) {
        lo[0] = 0.0; lo[1] = -all[1][0]; lo[2] = 0.0;
        lamda2x(lo,lo);
        hi[0] = 0.0; hi[1] = all[1][1]; hi[2] = 0.0;
        lamda2x(hi,hi);
        if (boundary[1][0] == 2) boxlo[1] = lo[1] - small[1];
        else if (boundary[1][0] == 3) boxlo[1] = MIN(lo[1]-small[1],minylo);
        if (boundary[1][1] == 2) boxhi[1] = hi[1] + small[1];
        else if (boundary[1][1] == 3) boxhi[1] = MAX(hi[1]+small[1],minyhi);
        if (boxlo[1] > boxhi[1]) error->all(FLERR,"Illegal simulation box");
      }
      if (zperiodic == 0) {
        lo[0] = 0.0; lo[1] = 0.0; lo[2] = -all[2][0];
        lamda2x(lo,lo);
        hi[0] = 0.0; hi[1] = 0.0; hi[2] = all[2][1];
        lamda2x(hi,hi);
        if (boundary[2][0] == 2) boxlo[2] = lo[2] - small[2];
        else if (boundary[2][0] == 3) boxlo[2] = MIN(lo[2]-small[2],minzlo);
        if (boundary[2][1] == 2) boxhi[2] = hi[2] + small[2];
        else if (boundary[2][1] == 3) boxhi[2] = MAX(hi[2]+small[2],minzhi);
        if (boxlo[2] > boxhi[2]) error->all(FLERR,"Illegal simulation box");
      }
    }
  }
  set_global_box();
  set_local_box();
  if (nonperiodic == 2 && triclinic) {
    x2lamda(atom->nlocal);
    pbc();
  }
}
void Domain::pbc()
{
  int nlocal = atom->nlocal;
  if (!nlocal) return;
  int i;
  imageint idim,otherdims;
  double *lo,*hi,*period;
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  imageint *image = atom->image;
  double *coord;
  int n3 = 3*nlocal;
  coord = &x[0][0];
  int flag = 0;
  for (i = 0; i < n3; i++)
    if (!std::isfinite(*coord++)) flag = 1;
  if (flag) error->one(FLERR,"Non-numeric atom coords - simulation unstable");
  if (triclinic == 0) {
    lo = boxlo;
    hi = boxhi;
    period = prd;
  } else {
    lo = boxlo_lamda;
    hi = boxhi_lamda;
    period = prd_lamda;
  }
  for (i = 0; i < nlocal; i++) {
    if (xperiodic) {
      if (x[i][0] < lo[0]) {
        x[i][0] += period[0];
        if (deform_vremap && mask[i] & deform_groupbit) v[i][0] += h_rate[0];
        idim = image[i] & IMGMASK;
        otherdims = image[i] ^ idim;
        idim--;
        idim &= IMGMASK;
        image[i] = otherdims | idim;
      }
      if (x[i][0] >= hi[0]) {
        x[i][0] -= period[0];
        x[i][0] = MAX(x[i][0],lo[0]);
        if (deform_vremap && mask[i] & deform_groupbit) v[i][0] -= h_rate[0];
        idim = image[i] & IMGMASK;
        otherdims = image[i] ^ idim;
        idim++;
        idim &= IMGMASK;
        image[i] = otherdims | idim;
      }
    }
    if (yperiodic) {
      if (x[i][1] < lo[1]) {
        x[i][1] += period[1];
        if (deform_vremap && mask[i] & deform_groupbit) {
          v[i][0] += h_rate[5];
          v[i][1] += h_rate[1];
        }
        idim = (image[i] >> IMGBITS) & IMGMASK;
        otherdims = image[i] ^ (idim << IMGBITS);
        idim--;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMGBITS);
      }
      if (x[i][1] >= hi[1]) {
        x[i][1] -= period[1];
        x[i][1] = MAX(x[i][1],lo[1]);
        if (deform_vremap && mask[i] & deform_groupbit) {
          v[i][0] -= h_rate[5];
          v[i][1] -= h_rate[1];
        }
        idim = (image[i] >> IMGBITS) & IMGMASK;
        otherdims = image[i] ^ (idim << IMGBITS);
        idim++;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMGBITS);
      }
    }
    if (zperiodic) {
      if (x[i][2] < lo[2]) {
        x[i][2] += period[2];
        if (deform_vremap && mask[i] & deform_groupbit) {
          v[i][0] += h_rate[4];
          v[i][1] += h_rate[3];
          v[i][2] += h_rate[2];
        }
        idim = image[i] >> IMG2BITS;
        otherdims = image[i] ^ (idim << IMG2BITS);
        idim--;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMG2BITS);
      }
      if (x[i][2] >= hi[2]) {
        x[i][2] -= period[2];
        x[i][2] = MAX(x[i][2],lo[2]);
        if (deform_vremap && mask[i] & deform_groupbit) {
          v[i][0] -= h_rate[4];
          v[i][1] -= h_rate[3];
          v[i][2] -= h_rate[2];
        }
        idim = image[i] >> IMG2BITS;
        otherdims = image[i] ^ (idim << IMG2BITS);
        idim++;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMG2BITS);
      }
    }
  }
}
int Domain::inside(double* x)
{
  double *lo,*hi;
  double lamda[3];
  if (triclinic == 0) {
    lo = boxlo;
    hi = boxhi;
    if (x[0] < lo[0] || x[0] >= hi[0] ||
        x[1] < lo[1] || x[1] >= hi[1] ||
        x[2] < lo[2] || x[2] >= hi[2]) return 0;
    else return 1;
  } else {
    lo = boxlo_lamda;
    hi = boxhi_lamda;
    x2lamda(x,lamda);
    if (lamda[0] < lo[0] || lamda[0] >= hi[0] ||
        lamda[1] < lo[1] || lamda[1] >= hi[1] ||
        lamda[2] < lo[2] || lamda[2] >= hi[2]) return 0;
    else return 1;
  }
}
int Domain::inside_nonperiodic(double* x)
{
  double *lo,*hi;
  double lamda[3];
  if (xperiodic && yperiodic && zperiodic) return 1;
  if (triclinic == 0) {
    lo = boxlo;
    hi = boxhi;
    if (!xperiodic && (x[0] < lo[0] || x[0] >= hi[0])) return 0;
    if (!yperiodic && (x[1] < lo[1] || x[1] >= hi[1])) return 0;
    if (!zperiodic && (x[2] < lo[2] || x[2] >= hi[2])) return 0;
    return 1;
  } else {
    lo = boxlo_lamda;
    hi = boxhi_lamda;
    x2lamda(x,lamda);
    if (!xperiodic && (lamda[0] < lo[0] || lamda[0] >= hi[0])) return 0;
    if (!yperiodic && (lamda[1] < lo[1] || lamda[1] >= hi[1])) return 0;
    if (!zperiodic && (lamda[2] < lo[2] || lamda[2] >= hi[2])) return 0;
    return 1;
  }
}
void Domain::subbox_too_small_check(double thresh)
{
  int flag = 0;
  if (!triclinic) {
    if (subhi[0]-sublo[0] < thresh || subhi[1]-sublo[1] < thresh) flag = 1;
    if (dimension == 3 && subhi[2]-sublo[2] < thresh) flag = 1;
  } else {
    double delta = subhi_lamda[0] - sublo_lamda[0];
    if (delta*prd[0] < thresh) flag = 1;
    delta = subhi_lamda[1] - sublo_lamda[1];
    if (delta*prd[1] < thresh) flag = 1;
    if (dimension == 3) {
      delta = subhi_lamda[2] - sublo_lamda[2];
      if (delta*prd[2] < thresh) flag = 1;
    }
  }
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall && comm->me == 0)
    error->warning(FLERR,"Proc sub-domain size < neighbor skin, "
                   "could lead to lost atoms");
}
static constexpr double MAXIMGCOUNT = 16;
void Domain::minimum_image(double &dx, double &dy, double &dz) const
{
  if (triclinic == 0) {
    if (xperiodic) {
      if (fabs(dx) > (MAXIMGCOUNT * xprd))
        error->one(FLERR, "Atoms have moved too far apart ({}) for minimum image\n", dx);
      while (fabs(dx) > xprd_half) {
        if (dx < 0.0) dx += xprd;
        else dx -= xprd;
      }
    }
    if (yperiodic) {
      if (fabs(dy) > (MAXIMGCOUNT * yprd))
        error->one(FLERR, "Atoms have moved too far apart ({}) for minimum image\n", dy);
      while (fabs(dy) > yprd_half) {
        if (dy < 0.0) dy += yprd;
        else dy -= yprd;
      }
    }
    if (zperiodic) {
      if (fabs(dz) > (MAXIMGCOUNT * zprd))
        error->one(FLERR, "Atoms have moved too far apart ({}) for minimum image\n", dz);
      while (fabs(dz) > zprd_half) {
        if (dz < 0.0) dz += zprd;
        else dz -= zprd;
      }
    }
  } else {
    if (zperiodic) {
      if (fabs(dz) > (MAXIMGCOUNT * zprd))
        error->one(FLERR, "Atoms have moved too far apart ({}) for minimum image\n", dz);
      while (fabs(dz) > zprd_half) {
        if (dz < 0.0) {
          dz += zprd;
          dy += yz;
          dx += xz;
        } else {
          dz -= zprd;
          dy -= yz;
          dx -= xz;
        }
      }
    }
    if (yperiodic) {
      if (fabs(dy) > (MAXIMGCOUNT * yprd))
        error->one(FLERR, "Atoms have moved too far apart ({}) for minimum image\n", dy);
      while (fabs(dy) > yprd_half) {
        if (dy < 0.0) {
          dy += yprd;
          dx += xy;
        } else {
          dy -= yprd;
          dx -= xy;
        }
      }
    }
    if (xperiodic) {
      if (fabs(dx) > (MAXIMGCOUNT * xprd))
        error->one(FLERR, "Atoms have moved too far apart ({}) for minimum image\n", dx);
      while (fabs(dx) > xprd_half) {
        if (dx < 0.0) dx += xprd;
        else dx -= xprd;
      }
    }
  }
}
int Domain::closest_image(int i, int j)
{
  if (j < 0) return j;
  int *sametag = atom->sametag;
  double **x = atom->x;
  double *xi = x[i];
  int closest = j;
  double delx = xi[0] - x[j][0];
  double dely = xi[1] - x[j][1];
  double delz = xi[2] - x[j][2];
  double rsqmin = delx*delx + dely*dely + delz*delz;
  double rsq;
  while (sametag[j] >= 0) {
    j = sametag[j];
    delx = xi[0] - x[j][0];
    dely = xi[1] - x[j][1];
    delz = xi[2] - x[j][2];
    rsq = delx*delx + dely*dely + delz*delz;
    if (rsq < rsqmin) {
      rsqmin = rsq;
      closest = j;
    }
  }
  return closest;
}
int Domain::closest_image(const double * const pos, int j)
{
  if (j < 0) return j;
  const int * const sametag = atom->sametag;
  const double * const * const x = atom->x;
  int closest = j;
  double delx = pos[0] - x[j][0];
  double dely = pos[1] - x[j][1];
  double delz = pos[2] - x[j][2];
  double rsqmin = delx*delx + dely*dely + delz*delz;
  double rsq;
  while (sametag[j] >= 0) {
    j = sametag[j];
    delx = pos[0] - x[j][0];
    dely = pos[1] - x[j][1];
    delz = pos[2] - x[j][2];
    rsq = delx*delx + dely*dely + delz*delz;
    if (rsq < rsqmin) {
      rsqmin = rsq;
      closest = j;
    }
  }
  return closest;
}
void Domain::closest_image(const double * const xi, const double * const xj, double * const xjimage)
{
  double dx = xj[0] - xi[0];
  double dy = xj[1] - xi[1];
  double dz = xj[2] - xi[2];
  if (triclinic == 0) {
    if (xperiodic) {
      if (dx < 0.0) {
        while (dx < 0.0) dx += xprd;
        if (dx > xprd_half) dx -= xprd;
      } else {
        while (dx > 0.0) dx -= xprd;
        if (dx < -xprd_half) dx += xprd;
      }
    }
    if (yperiodic) {
      if (dy < 0.0) {
        while (dy < 0.0) dy += yprd;
        if (dy > yprd_half) dy -= yprd;
      } else {
        while (dy > 0.0) dy -= yprd;
        if (dy < -yprd_half) dy += yprd;
      }
    }
    if (zperiodic) {
      if (dz < 0.0) {
        while (dz < 0.0) dz += zprd;
        if (dz > zprd_half) dz -= zprd;
      } else {
        while (dz > 0.0) dz -= zprd;
        if (dz < -zprd_half) dz += zprd;
      }
    }
  } else {
    if (zperiodic) {
      if (dz < 0.0) {
        while (dz < 0.0) {
          dz += zprd;
          dy += yz;
          dx += xz;
        }
        if (dz > zprd_half) {
          dz -= zprd;
          dy -= yz;
          dx -= xz;
        }
      } else {
        while (dz > 0.0) {
          dz -= zprd;
          dy -= yz;
          dx -= xz;
        }
        if (dz < -zprd_half) {
          dz += zprd;
          dy += yz;
          dx += xz;
        }
      }
    }
    if (yperiodic) {
      if (dy < 0.0) {
        while (dy < 0.0) {
          dy += yprd;
          dx += xy;
        }
        if (dy > yprd_half) {
          dy -= yprd;
          dx -= xy;
        }
      } else {
        while (dy > 0.0) {
          dy -= yprd;
          dx -= xy;
        }
        if (dy < -yprd_half) {
          dy += yprd;
          dx += xy;
        }
      }
    }
    if (xperiodic) {
      if (dx < 0.0) {
        while (dx < 0.0) dx += xprd;
        if (dx > xprd_half) dx -= xprd;
      } else {
        while (dx > 0.0) dx -= xprd;
        if (dx < -xprd_half) dx += xprd;
      }
    }
  }
  xjimage[0] = xi[0] + dx;
  xjimage[1] = xi[1] + dy;
  xjimage[2] = xi[2] + dz;
}
void Domain::remap(double *x, imageint &image)
{
  double *lo,*hi,*period,*coord;
  double lamda[3];
  imageint idim,otherdims;
  if (triclinic == 0) {
    lo = boxlo;
    hi = boxhi;
    period = prd;
    coord = x;
  } else {
    lo = boxlo_lamda;
    hi = boxhi_lamda;
    period = prd_lamda;
    x2lamda(x,lamda);
    coord = lamda;
  }
  if (xperiodic) {
    while (coord[0] < lo[0]) {
      coord[0] += period[0];
      idim = image & IMGMASK;
      otherdims = image ^ idim;
      idim--;
      idim &= IMGMASK;
      image = otherdims | idim;
    }
    while (coord[0] >= hi[0]) {
      coord[0] -= period[0];
      idim = image & IMGMASK;
      otherdims = image ^ idim;
      idim++;
      idim &= IMGMASK;
      image = otherdims | idim;
    }
    coord[0] = MAX(coord[0],lo[0]);
  }
  if (yperiodic) {
    while (coord[1] < lo[1]) {
      coord[1] += period[1];
      idim = (image >> IMGBITS) & IMGMASK;
      otherdims = image ^ (idim << IMGBITS);
      idim--;
      idim &= IMGMASK;
      image = otherdims | (idim << IMGBITS);
    }
    while (coord[1] >= hi[1]) {
      coord[1] -= period[1];
      idim = (image >> IMGBITS) & IMGMASK;
      otherdims = image ^ (idim << IMGBITS);
      idim++;
      idim &= IMGMASK;
      image = otherdims | (idim << IMGBITS);
    }
    coord[1] = MAX(coord[1],lo[1]);
  }
  if (zperiodic) {
    while (coord[2] < lo[2]) {
      coord[2] += period[2];
      idim = image >> IMG2BITS;
      otherdims = image ^ (idim << IMG2BITS);
      idim--;
      idim &= IMGMASK;
      image = otherdims | (idim << IMG2BITS);
    }
    while (coord[2] >= hi[2]) {
      coord[2] -= period[2];
      idim = image >> IMG2BITS;
      otherdims = image ^ (idim << IMG2BITS);
      idim++;
      idim &= IMGMASK;
      image = otherdims | (idim << IMG2BITS);
    }
    coord[2] = MAX(coord[2],lo[2]);
  }
  if (triclinic) lamda2x(coord,x);
}
void Domain::remap(double *x)
{
  double *lo,*hi,*period,*coord;
  double lamda[3];
  if (triclinic == 0) {
    lo = boxlo;
    hi = boxhi;
    period = prd;
    coord = x;
  } else {
    lo = boxlo_lamda;
    hi = boxhi_lamda;
    period = prd_lamda;
    x2lamda(x,lamda);
    coord = lamda;
  }
  if (xperiodic) {
    while (coord[0] < lo[0]) coord[0] += period[0];
    while (coord[0] >= hi[0]) coord[0] -= period[0];
    coord[0] = MAX(coord[0],lo[0]);
  }
  if (yperiodic) {
    while (coord[1] < lo[1]) coord[1] += period[1];
    while (coord[1] >= hi[1]) coord[1] -= period[1];
    coord[1] = MAX(coord[1],lo[1]);
  }
  if (zperiodic) {
    while (coord[2] < lo[2]) coord[2] += period[2];
    while (coord[2] >= hi[2]) coord[2] -= period[2];
    coord[2] = MAX(coord[2],lo[2]);
  }
  if (triclinic) lamda2x(coord,x);
}
void Domain::remap_near(double *xnew, double *xold)
{
  int n;
  double *coordnew,*coordold,*period,*half;
  double lamdanew[3],lamdaold[3];
  if (triclinic == 0) {
    period = prd;
    half = prd_half;
    coordnew = xnew;
    coordold = xold;
  } else {
    period = prd_lamda;
    half = prd_half_lamda;
    x2lamda(xnew,lamdanew);
    coordnew = lamdanew;
    x2lamda(xold,lamdaold);
    coordold = lamdaold;
  }
  if (xperiodic) {
    if (coordnew[0]-coordold[0] > period[0]) {
      n = static_cast<int> ((coordnew[0]-coordold[0])/period[0]);
      coordnew[0] -= n*period[0];
    }
    while (coordnew[0]-coordold[0] > half[0]) coordnew[0] -= period[0];
    if (coordold[0]-coordnew[0] > period[0]) {
      n = static_cast<int> ((coordold[0]-coordnew[0])/period[0]);
      coordnew[0] += n*period[0];
    }
    while (coordold[0]-coordnew[0] > half[0]) coordnew[0] += period[0];
  }
  if (yperiodic) {
    if (coordnew[1]-coordold[1] > period[1]) {
      n = static_cast<int> ((coordnew[1]-coordold[1])/period[1]);
      coordnew[1] -= n*period[1];
    }
    while (coordnew[1]-coordold[1] > half[1]) coordnew[1] -= period[1];
    if (coordold[1]-coordnew[1] > period[1]) {
      n = static_cast<int> ((coordold[1]-coordnew[1])/period[1]);
      coordnew[1] += n*period[1];
    }
    while (coordold[1]-coordnew[1] > half[1]) coordnew[1] += period[1];
  }
  if (zperiodic) {
    if (coordnew[2]-coordold[2] > period[2]) {
      n = static_cast<int> ((coordnew[2]-coordold[2])/period[2]);
      coordnew[2] -= n*period[2];
    }
    while (coordnew[2]-coordold[2] > half[2]) coordnew[2] -= period[2];
    if (coordold[2]-coordnew[2] > period[2]) {
      n = static_cast<int> ((coordold[2]-coordnew[2])/period[2]);
      coordnew[2] += n*period[2];
    }
    while (coordold[2]-coordnew[2] > half[2]) coordnew[2] += period[2];
  }
  if (triclinic) lamda2x(coordnew,xnew);
}
void Domain::unmap(double *x, imageint image)
{
  int xbox = (image & IMGMASK) - IMGMAX;
  int ybox = (image >> IMGBITS & IMGMASK) - IMGMAX;
  int zbox = (image >> IMG2BITS) - IMGMAX;
  if (triclinic == 0) {
    x[0] += xbox*xprd;
    x[1] += ybox*yprd;
    x[2] += zbox*zprd;
  } else {
    x[0] += h[0]*xbox + h[5]*ybox + h[4]*zbox;
    x[1] += h[1]*ybox + h[3]*zbox;
    x[2] += h[2]*zbox;
  }
}
void Domain::unmap(const double *x, imageint image, double *y)
{
  int xbox = (image & IMGMASK) - IMGMAX;
  int ybox = (image >> IMGBITS & IMGMASK) - IMGMAX;
  int zbox = (image >> IMG2BITS) - IMGMAX;
  if (triclinic == 0) {
    y[0] = x[0] + xbox*xprd;
    y[1] = x[1] + ybox*yprd;
    y[2] = x[2] + zbox*zprd;
  } else {
    y[0] = x[0] + h[0]*xbox + h[5]*ybox + h[4]*zbox;
    y[1] = x[1] + h[1]*ybox + h[3]*zbox;
    y[2] = x[2] + h[2]*zbox;
  }
}
void Domain::image_flip(int m, int n, int p)
{
  imageint *image = atom->image;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    int xbox = (image[i] & IMGMASK) - IMGMAX;
    int ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
    int zbox = (image[i] >> IMG2BITS) - IMGMAX;
    ybox -= p*zbox;
    xbox -= m*ybox + n*zbox;
    image[i] = ((imageint) (xbox + IMGMAX) & IMGMASK) |
      (((imageint) (ybox + IMGMAX) & IMGMASK) << IMGBITS) |
      (((imageint) (zbox + IMGMAX) & IMGMASK) << IMG2BITS);
  }
}
int Domain::ownatom(int , double *x, imageint *image, int shrinkexceed)
{
  double lamda[3];
  double *coord,*blo,*bhi,*slo,*shi;
  if (image) remap(x,*image);
  else remap(x);
  if (triclinic) {
    x2lamda(x,lamda);
    if (xperiodic && (lamda[0] < 0.0 || lamda[0] >= 1.0)) lamda[0] = 0.0;
    if (yperiodic && (lamda[1] < 0.0 || lamda[1] >= 1.0)) lamda[1] = 0.0;
    if (zperiodic && (lamda[2] < 0.0 || lamda[2] >= 1.0)) lamda[2] = 0.0;
    coord = lamda;
  } else coord = x;
  if (triclinic == 0) {
    blo = boxlo;
    bhi = boxhi;
    slo = sublo;
    shi = subhi;
  } else {
    blo = boxlo_lamda;
    bhi = boxhi_lamda;
    slo = sublo_lamda;
    shi = subhi_lamda;
  }
  if (coord[0] >= slo[0] && coord[0] < shi[0] &&
      coord[1] >= slo[1] && coord[1] < shi[1] &&
      coord[2] >= slo[2] && coord[2] < shi[2]) return 1;
  if (shrinkexceed) {
    int outside = 0;
    if (coord[0] < blo[0] && boundary[0][0] > 1) outside = 1;
    if (coord[0] >= bhi[0] && boundary[0][1] > 1) outside = 1;
    if (coord[1] < blo[1] && boundary[1][0] > 1) outside = 1;
    if (coord[1] >= bhi[1] && boundary[1][1] > 1) outside = 1;
    if (coord[2] < blo[2] && boundary[2][0] > 1) outside = 1;
    if (coord[2] >= bhi[2] && boundary[2][1] > 1) outside = 1;
    if (!outside) return 0;
    double newcoord[3];
    if (coord[0] < blo[0] && boundary[0][0] > 1) newcoord[0] = blo[0];
    else if (coord[0] >= bhi[0] && boundary[0][1] > 1) newcoord[0] = bhi[0];
    else newcoord[0] = coord[0];
    if (coord[1] < blo[1] && boundary[1][0] > 1) newcoord[1] = blo[1];
    else if (coord[1] >= bhi[1] && boundary[1][1] > 1) newcoord[1] = bhi[1];
    else newcoord[1] = coord[1];
    if (coord[2] < blo[2] && boundary[2][0] > 1) newcoord[2] = blo[2];
    else if (coord[2] >= bhi[2] && boundary[2][1] > 1) newcoord[2] = bhi[2];
    else newcoord[2] = coord[2];
    if (newcoord[0] >= slo[0] && newcoord[0] <= shi[0] &&
        newcoord[1] >= slo[1] && newcoord[1] <= shi[1] &&
        newcoord[2] >= slo[2] && newcoord[2] <= shi[2]) return 1;
  }
  return 0;
}
void Domain::set_lattice(int narg, char **arg)
{
  if (lattice) delete lattice;
  lattice = nullptr;
  lattice = new Lattice(lmp,narg,arg);
}
void Domain::add_region(int narg, char **arg)
{
  if (narg < 2) utils::missing_cmd_args(FLERR, "region", error);
  if (strcmp(arg[1],"none") == 0)
    error->all(FLERR,"Unrecognized region style 'none'");
  if (get_region_by_id(arg[0])) error->all(FLERR,"Reuse of region ID {}", arg[0]);
  Region *newregion = nullptr;
  if (!newregion && (region_map->find(arg[1]) != region_map->end())) {
    RegionCreator &region_creator = (*region_map)[arg[1]];
    newregion = region_creator(lmp, narg, arg);
  }
  newregion->init();
  regions.insert(newregion);
}
Region *Domain::get_region_by_id(const std::string &name) const
{
  for (auto &reg : regions)
    if (name == reg->id) return reg;
  return nullptr;
}
const std::vector<Region *> Domain::get_region_by_style(const std::string &name) const
{
  std::vector<Region *> matches;
  if (name.empty()) return matches;
  for (auto &reg : regions)
    if (name == reg->style) matches.push_back(reg);
  return matches;
}
const std::vector<Region *> Domain::get_region_list()
{
  return std::vector<Region *>(regions.begin(), regions.end());
}
void Domain::set_boundary(int narg, char **arg, int flag)
{
  if (narg != 3) error->all(FLERR,"Illegal boundary command: expected 3 arguments but found {}", narg);
  char c;
  for (int idim = 0; idim < 3; idim++)
    for (int iside = 0; iside < 2; iside++) {
      if (iside == 0) c = arg[idim][0];
      else if (iside == 1 && strlen(arg[idim]) == 1) c = arg[idim][0];
      else c = arg[idim][1];
      if (c == 'p') boundary[idim][iside] = 0;
      else if (c == 'f') boundary[idim][iside] = 1;
      else if (c == 's') boundary[idim][iside] = 2;
      else if (c == 'm') boundary[idim][iside] = 3;
      else {
        if (flag == 0) error->all(FLERR,"Unknown boundary keyword: {}", c);
        if (flag == 1) error->all(FLERR,"Unknown change_box keyword: {}", c);
      }
    }
  for (int idim = 0; idim < 3; idim++)
    if ((boundary[idim][0] == 0 && boundary[idim][1]) ||
        (boundary[idim][0] && boundary[idim][1] == 0))
      error->all(FLERR,"Both sides of boundary must be periodic");
  if (boundary[0][0] == 0) xperiodic = 1;
  else xperiodic = 0;
  if (boundary[1][0] == 0) yperiodic = 1;
  else yperiodic = 0;
  if (boundary[2][0] == 0) zperiodic = 1;
  else zperiodic = 0;
  int pflag=0;
  if ((periodicity[0] && !xperiodic)
      || (periodicity[1] && !yperiodic)
      || (periodicity[2] && !zperiodic)) pflag = 1;
  periodicity[0] = xperiodic;
  periodicity[1] = yperiodic;
  periodicity[2] = zperiodic;
  nonperiodic = 0;
  if (xperiodic == 0 || yperiodic == 0 || zperiodic == 0) {
    nonperiodic = 1;
    if (boundary[0][0] >= 2 || boundary[0][1] >= 2 ||
        boundary[1][0] >= 2 || boundary[1][1] >= 2 ||
        boundary[2][0] >= 2 || boundary[2][1] >= 2) nonperiodic = 2;
  }
  if (pflag) {
    pflag = 0;
    for (int i=0; i < atom->nlocal; ++i) {
      int xbox = (atom->image[i] & IMGMASK) - IMGMAX;
      int ybox = (atom->image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      int zbox = (atom->image[i] >> IMG2BITS) - IMGMAX;
      if ((!xperiodic) && (xbox != 0)) { xbox = 0; pflag = 1; }
      if ((!yperiodic) && (ybox != 0)) { ybox = 0; pflag = 1; }
      if ((!zperiodic) && (zbox != 0)) { zbox = 0; pflag = 1; }
      atom->image[i] = ((imageint) (xbox + IMGMAX) & IMGMASK) |
        (((imageint) (ybox + IMGMAX) & IMGMASK) << IMGBITS) |
        (((imageint) (zbox + IMGMAX) & IMGMASK) << IMG2BITS);
    }
    int flag_all;
    MPI_Allreduce(&pflag,&flag_all, 1, MPI_INT, MPI_SUM, world);
    if ((flag_all > 0) && (comm->me == 0))
      error->warning(FLERR,"Resetting image flags for non-periodic dimensions");
  }
}
void Domain::print_box(const std::string &prefix)
{
  if (comm->me == 0) {
    std::string mesg = prefix;
    if (triclinic == 0) {
      mesg += fmt::format("orthogonal box = ({:.8} {:.8} {:.8}) to "
                          "({:.8} {:.8} {:.8})\n",boxlo[0],boxlo[1],
                          boxlo[2],boxhi[0],boxhi[1],boxhi[2]);
    } else {
      mesg += fmt::format("triclinic box = ({:.8} {:.8} {:.8}) to "
                          "({:.8} {:.8} {:.8}) with tilt "
                          "({:.8} {:.8} {:.8})\n",boxlo[0],boxlo[1],
                          boxlo[2],boxhi[0],boxhi[1],boxhi[2],xy,xz,yz);
    }
    utils::logmesg(lmp,mesg);
  }
}
void Domain::boundary_string(char *str)
{
  int m = 0;
  for (int idim = 0; idim < 3; idim++) {
    for (int iside = 0; iside < 2; iside++) {
      if (boundary[idim][iside] == 0) str[m++] = 'p';
      else if (boundary[idim][iside] == 1) str[m++] = 'f';
      else if (boundary[idim][iside] == 2) str[m++] = 's';
      else if (boundary[idim][iside] == 3) str[m++] = 'm';
    }
    str[m++] = ' ';
  }
  str[8] = '\0';
}
void Domain::lamda2x(int n)
{
  double **x = atom->x;
  for (int i = 0; i < n; i++) {
    x[i][0] = h[0]*x[i][0] + h[5]*x[i][1] + h[4]*x[i][2] + boxlo[0];
    x[i][1] = h[1]*x[i][1] + h[3]*x[i][2] + boxlo[1];
    x[i][2] = h[2]*x[i][2] + boxlo[2];
  }
}
void Domain::x2lamda(int n)
{
  double delta[3];
  double **x = atom->x;
  for (int i = 0; i < n; i++) {
    delta[0] = x[i][0] - boxlo[0];
    delta[1] = x[i][1] - boxlo[1];
    delta[2] = x[i][2] - boxlo[2];
    x[i][0] = h_inv[0]*delta[0] + h_inv[5]*delta[1] + h_inv[4]*delta[2];
    x[i][1] = h_inv[1]*delta[1] + h_inv[3]*delta[2];
    x[i][2] = h_inv[2]*delta[2];
  }
}
void Domain::lamda2x(double *lamda, double *x)
{
  x[0] = h[0]*lamda[0] + h[5]*lamda[1] + h[4]*lamda[2] + boxlo[0];
  x[1] = h[1]*lamda[1] + h[3]*lamda[2] + boxlo[1];
  x[2] = h[2]*lamda[2] + boxlo[2];
}
void Domain::x2lamda(double *x, double *lamda)
{
  double delta[3];
  delta[0] = x[0] - boxlo[0];
  delta[1] = x[1] - boxlo[1];
  delta[2] = x[2] - boxlo[2];
  lamda[0] = h_inv[0]*delta[0] + h_inv[5]*delta[1] + h_inv[4]*delta[2];
  lamda[1] = h_inv[1]*delta[1] + h_inv[3]*delta[2];
  lamda[2] = h_inv[2]*delta[2];
}
void Domain::x2lamda(double *x, double *lamda,
                     double *my_boxlo, double *my_h_inv)
{
  double delta[3];
  delta[0] = x[0] - my_boxlo[0];
  delta[1] = x[1] - my_boxlo[1];
  delta[2] = x[2] - my_boxlo[2];
  lamda[0] = my_h_inv[0]*delta[0] + my_h_inv[5]*delta[1] + my_h_inv[4]*delta[2];
  lamda[1] = my_h_inv[1]*delta[1] + my_h_inv[3]*delta[2];
  lamda[2] = my_h_inv[2]*delta[2];
}
void Domain::bbox(double *lo, double *hi, double *bboxlo, double *bboxhi)
{
  double x[3];
  bboxlo[0] = bboxlo[1] = bboxlo[2] = BIG;
  bboxhi[0] = bboxhi[1] = bboxhi[2] = -BIG;
  x[0] = lo[0]; x[1] = lo[1]; x[2] = lo[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);
  x[0] = hi[0]; x[1] = lo[1]; x[2] = lo[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);
  x[0] = lo[0]; x[1] = hi[1]; x[2] = lo[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);
  x[0] = hi[0]; x[1] = hi[1]; x[2] = lo[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);
  x[0] = lo[0]; x[1] = lo[1]; x[2] = hi[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);
  x[0] = hi[0]; x[1] = lo[1]; x[2] = hi[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);
  x[0] = lo[0]; x[1] = hi[1]; x[2] = hi[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);
  x[0] = hi[0]; x[1] = hi[1]; x[2] = hi[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);
}
void Domain::box_corners()
{
  lamda_box_corners(boxlo_lamda,boxhi_lamda);
}
void Domain::subbox_corners()
{
  lamda_box_corners(sublo_lamda,subhi_lamda);
}
void Domain::lamda_box_corners(double *lo, double *hi)
{
  corners[0][0] = lo[0]; corners[0][1] = lo[1]; corners[0][2] = lo[2];
  lamda2x(corners[0],corners[0]);
  corners[1][0] = hi[0]; corners[1][1] = lo[1]; corners[1][2] = lo[2];
  lamda2x(corners[1],corners[1]);
  corners[2][0] = lo[0]; corners[2][1] = hi[1]; corners[2][2] = lo[2];
  lamda2x(corners[2],corners[2]);
  corners[3][0] = hi[0]; corners[3][1] = hi[1]; corners[3][2] = lo[2];
  lamda2x(corners[3],corners[3]);
  corners[4][0] = lo[0]; corners[4][1] = lo[1]; corners[4][2] = hi[2];
  lamda2x(corners[4],corners[4]);
  corners[5][0] = hi[0]; corners[5][1] = lo[1]; corners[5][2] = hi[2];
  lamda2x(corners[5],corners[5]);
  corners[6][0] = lo[0]; corners[6][1] = hi[1]; corners[6][2] = hi[2];
  lamda2x(corners[6],corners[6]);
  corners[7][0] = hi[0]; corners[7][1] = hi[1]; corners[7][2] = hi[2];
  lamda2x(corners[7],corners[7]);
}
