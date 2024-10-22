#include "domain.h"
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
#include "style_region.h"
#include "universe.h"
#include "update.h"
#include <cmath>
#include <cstring>
using namespace LAMMPS_NS;
#define BIG 1.0e20
#define SMALL 1.0e-4
#define DELTAREGION 4
#define BONDSTRETCH 1.1
template <typename T>
static Region *region_creator(LAMMPS *lmp, int narg, char **arg) {
  return new T(lmp, narg, arg);
}
Domain::Domain(LAMMPS *lmp) : Pointers(lmp) {
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
  h_rate[0] = h_rate[1] = h_rate[2] = h_rate[3] = h_rate[4] = h_rate[5] = 0.0;
  h_ratelo[0] = h_ratelo[1] = h_ratelo[2] = 0.0;
  prd_lamda[0] = prd_lamda[1] = prd_lamda[2] = 1.0;
  prd_half_lamda[0] = prd_half_lamda[1] = prd_half_lamda[2] = 0.5;
  boxlo_lamda[0] = boxlo_lamda[1] = boxlo_lamda[2] = 0.0;
  boxhi_lamda[0] = boxhi_lamda[1] = boxhi_lamda[2] = 1.0;
  lattice = nullptr;
  auto args = new char *[2];
  args[0] = (char *)"none";
  args[1] = (char *)"1.0";
  set_lattice(2, args);
  delete[] args;
  copymode = 0;
  region_map = new RegionCreatorMap();
#define REGION_CLASS
#define RegionStyle(key, Class) (*region_map)[#key] = &region_creator<Class>;
#include "style_region.h"
#undef RegionStyle
#undef REGION_CLASS
}
Domain::~Domain() {
  if (copymode)
    return;
  for (auto &reg : regions)
    delete reg;
  regions.clear();
  delete lattice;
  delete region_map;
}
void Domain::init() {
  box_change_size = box_change_shape = box_change_domain = 0;
  int box_change_x = 0, box_change_y = 0, box_change_z = 0;
  int box_change_yz = 0, box_change_xz = 0, box_change_xy = 0;
  const auto &fixes = modify->get_fix_list();
  if (nonperiodic == 2)
    box_change_size = 1;
  for (const auto &fix : fixes) {
    if (fix->box_change & Fix::BOX_CHANGE_SIZE)
      box_change_size = 1;
    if (fix->box_change & Fix::BOX_CHANGE_SHAPE)
      box_change_shape = 1;
    if (fix->box_change & Fix::BOX_CHANGE_DOMAIN)
      box_change_domain = 1;
    if (fix->box_change & Fix::BOX_CHANGE_X)
      box_change_x++;
    if (fix->box_change & Fix::BOX_CHANGE_Y)
      box_change_y++;
    if (fix->box_change & Fix::BOX_CHANGE_Z)
      box_change_z++;
    if (fix->box_change & Fix::BOX_CHANGE_YZ)
      box_change_yz++;
    if (fix->box_change & Fix::BOX_CHANGE_XZ)
      box_change_xz++;
    if (fix->box_change & Fix::BOX_CHANGE_XY)
      box_change_xy++;
  }
  std::string mesg = "Must not have multiple fixes change box parameter ";
#define CHECK_BOX_FIX_ERROR(par)                                               \
  if (box_change_##par > 1)                                                    \
  error->all(FLERR, (mesg + #par))
  CHECK_BOX_FIX_ERROR(x);
  CHECK_BOX_FIX_ERROR(y);
  CHECK_BOX_FIX_ERROR(z);
  CHECK_BOX_FIX_ERROR(yz);
  CHECK_BOX_FIX_ERROR(xz);
  CHECK_BOX_FIX_ERROR(xy);
#undef CHECK_BOX_FIX_ERROR
  box_change = 0;
  if (box_change_size || box_change_shape || box_change_domain)
    box_change = 1;
  deform_flag = deform_vremap = deform_groupbit = 0;
  for (auto &reg : regions)
    reg->init();
}
void Domain::set_initial_box(int expandflag) {
  small[0] = SMALL * (boxhi[0] - boxlo[0]);
  small[1] = SMALL * (boxhi[1] - boxlo[1]);
  small[2] = SMALL * (boxhi[2] - boxlo[2]);
  if (!expandflag)
    return;
  if (boundary[0][0] == 2)
    boxlo[0] -= small[0];
  else if (boundary[0][0] == 3)
    minxlo = boxlo[0];
  if (boundary[0][1] == 2)
    boxhi[0] += small[0];
  else if (boundary[0][1] == 3)
    minxhi = boxhi[0];
  if (boundary[1][0] == 2)
    boxlo[1] -= small[1];
  else if (boundary[1][0] == 3)
    minylo = boxlo[1];
  if (boundary[1][1] == 2)
    boxhi[1] += small[1];
  else if (boundary[1][1] == 3)
    minyhi = boxhi[1];
  if (boundary[2][0] == 2)
    boxlo[2] -= small[2];
  else if (boundary[2][0] == 3)
    minzlo = boxlo[2];
  if (boundary[2][1] == 2)
    boxhi[2] += small[2];
  else if (boundary[2][1] == 3)
    minzhi = boxhi[2];
}
void Domain::set_global_box() {
  prd[0] = xprd = boxhi[0] - boxlo[0];
  prd[1] = yprd = boxhi[1] - boxlo[1];
  prd[2] = zprd = boxhi[2] - boxlo[2];
  h[0] = xprd;
  h[1] = yprd;
  h[2] = zprd;
  h_inv[0] = 1.0 / h[0];
  h_inv[1] = 1.0 / h[1];
  h_inv[2] = 1.0 / h[2];
  prd_half[0] = xprd_half = 0.5 * xprd;
  prd_half[1] = yprd_half = 0.5 * yprd;
  prd_half[2] = zprd_half = 0.5 * zprd;
}
void Domain::set_lamda_box() {
  int *myloc = comm->myloc;
  double *xsplit = comm->xsplit;
  double *ysplit = comm->ysplit;
  double *zsplit = comm->zsplit;
  sublo_lamda[0] = xsplit[myloc[0]];
  subhi_lamda[0] = xsplit[myloc[0] + 1];
  sublo_lamda[1] = ysplit[myloc[1]];
  subhi_lamda[1] = ysplit[myloc[1] + 1];
  sublo_lamda[2] = zsplit[myloc[2]];
  subhi_lamda[2] = zsplit[myloc[2] + 1];
}
void Domain::set_local_box() {
  int *myloc = comm->myloc;
  int *procgrid = comm->procgrid;
  double *xsplit = comm->xsplit;
  double *ysplit = comm->ysplit;
  double *zsplit = comm->zsplit;
  sublo[0] = boxlo[0] + xprd * xsplit[myloc[0]];
  if (myloc[0] < procgrid[0] - 1)
    subhi[0] = boxlo[0] + xprd * xsplit[myloc[0] + 1];
  else
    subhi[0] = boxhi[0];
  sublo[1] = boxlo[1] + yprd * ysplit[myloc[1]];
  if (myloc[1] < procgrid[1] - 1)
    subhi[1] = boxlo[1] + yprd * ysplit[myloc[1] + 1];
  else
    subhi[1] = boxhi[1];
  sublo[2] = boxlo[2] + zprd * zsplit[myloc[2]];
  if (myloc[2] < procgrid[2] - 1)
    subhi[2] = boxlo[2] + zprd * zsplit[myloc[2] + 1];
  else
    subhi[2] = boxhi[2];
}
void Domain::reset_box() {
  if (atom->natoms == 0)
    return;
  set_global_box();
  set_local_box();
}
void Domain::pbc() {
  int nlocal = atom->nlocal;
  if (!nlocal)
    return;
  int i;
  imageint idim, otherdims;
  double *lo, *hi, *period;
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  imageint *image = atom->image;
  double *coord;
  int n3 = 3 * nlocal;
  coord = &x[0][0];
  int flag = 0;
  for (i = 0; i < n3; i++)
    if (!std::isfinite(*coord++))
      flag = 1;
  if (flag)
    error->one(FLERR, "Non-numeric atom coords - simulation unstable");
  lo = boxlo;
  hi = boxhi;
  period = prd;
  for (i = 0; i < nlocal; i++) {
    if (xperiodic) {
      if (x[i][0] < lo[0]) {
        x[i][0] += period[0];
        if (deform_vremap && mask[i] & deform_groupbit)
          v[i][0] += h_rate[0];
        idim = image[i] & IMGMASK;
        otherdims = image[i] ^ idim;
        idim--;
        idim &= IMGMASK;
        image[i] = otherdims | idim;
      }
      if (x[i][0] >= hi[0]) {
        x[i][0] -= period[0];
        x[i][0] = MAX(x[i][0], lo[0]);
        if (deform_vremap && mask[i] & deform_groupbit)
          v[i][0] -= h_rate[0];
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
        x[i][1] = MAX(x[i][1], lo[1]);
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
        x[i][2] = MAX(x[i][2], lo[2]);
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
void Domain::subbox_too_small_check(double thresh) {
  int flag = 0;
  if (!triclinic) {
    if (subhi[0] - sublo[0] < thresh || subhi[1] - sublo[1] < thresh)
      flag = 1;
    if (dimension == 3 && subhi[2] - sublo[2] < thresh)
      flag = 1;
  } else {
    double delta = subhi_lamda[0] - sublo_lamda[0];
    if (delta * prd[0] < thresh)
      flag = 1;
    delta = subhi_lamda[1] - sublo_lamda[1];
    if (delta * prd[1] < thresh)
      flag = 1;
    if (dimension == 3) {
      delta = subhi_lamda[2] - sublo_lamda[2];
      if (delta * prd[2] < thresh)
        flag = 1;
    }
  }
  int flagall;
  MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_SUM, world);
  if (flagall && comm->me == 0)
    error->warning(FLERR, "Proc sub-domain size < neighbor skin, "
                          "could lead to lost atoms");
}
static constexpr double MAXIMGCOUNT = 16;
void Domain::set_lattice(int narg, char **arg) {
  if (lattice)
    delete lattice;
  lattice = nullptr;
  lattice = new Lattice(lmp, narg, arg);
}
void Domain::add_region(int narg, char **arg) {
  if (narg < 2)
    utils::missing_cmd_args(FLERR, "region", error);
  if (strcmp(arg[1], "none") == 0)
    error->all(FLERR, "Unrecognized region style 'none'");
  if (get_region_by_id(arg[0]))
    error->all(FLERR, "Reuse of region ID {}", arg[0]);
  Region *newregion = nullptr;
  if (!newregion && (region_map->find(arg[1]) != region_map->end())) {
    RegionCreator &region_creator = (*region_map)[arg[1]];
    newregion = region_creator(lmp, narg, arg);
  }
  newregion->init();
  regions.insert(newregion);
}
Region *Domain::get_region_by_id(const std::string &name) const {
  for (auto &reg : regions)
    if (name == reg->id)
      return reg;
  return nullptr;
}
