#include "region.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "lattice.h"
#include "math_extra.h"
#include "update.h"
#include "variable.h"
#include <cmath>
#include <cstring>
using namespace LAMMPS_NS;
Region::Region(LAMMPS *lmp, int , char **arg) :
    Pointers(lmp), id(nullptr), style(nullptr), reglist(nullptr), contact(nullptr), xstr(nullptr),
    ystr(nullptr), zstr(nullptr), tstr(nullptr)
{
  id = utils::strdup(arg[0]);
  style = utils::strdup(arg[1]);
  varshape = 0;
  xstr = ystr = zstr = tstr = nullptr;
  dx = dy = dz = 0.0;
  size_restart = 5;
  Region::reset_vel();
  copymode = 0;
  nregion = 1;
}
Region::~Region()
{
  if (copymode) return;
  delete[] id;
  delete[] style;
  delete[] xstr;
  delete[] ystr;
  delete[] zstr;
  delete[] tstr;
}
void Region::init()
{
  vel_timestep = -1;
}
int Region::dynamic_check()
{
  if (dynamic || varshape) return 1;
  return 0;
}
void Region::prematch()
{
  if (varshape) shape_update();
  if (dynamic) pretransform();
}
int Region::match(double x, double y, double z)
{
  if (dynamic) inverse_transform(x, y, z);
  if (openflag) return 1;
  return !(inside(x, y, z) ^ interior);
}
void Region::add_contact(int n, double *x, double xp, double yp, double zp)
{
  double delx = x[0] - xp;
  double dely = x[1] - yp;
  double delz = x[2] - zp;
  contact[n].r = sqrt(delx * delx + dely * dely + delz * delz);
  contact[n].radius = 0;
  contact[n].delx = delx;
  contact[n].dely = dely;
  contact[n].delz = delz;
}
void Region::pretransform()
{
  if (moveflag) {
    if (xstr) dx = input->variable->compute_equal(xvar);
    if (ystr) dy = input->variable->compute_equal(yvar);
    if (zstr) dz = input->variable->compute_equal(zvar);
  }
  if (rotateflag) theta = input->variable->compute_equal(tvar);
}
void Region::forward_transform(double &x, double &y, double &z)
{
  if (rotateflag) rotate(x, y, z, theta);
  if (moveflag) {
    x += dx;
    y += dy;
    z += dz;
  }
}
void Region::inverse_transform(double &x, double &y, double &z)
{
  if (moveflag) {
    x -= dx;
    y -= dy;
    z -= dz;
  }
  if (rotateflag) rotate(x, y, z, -theta);
}
void Region::rotate(double &x, double &y, double &z, double angle)
{
  double a[3], b[3], c[3], d[3], disp[3];
  double sine = sin(angle);
  double cosine = cos(angle);
  d[0] = x - point[0];
  d[1] = y - point[1];
  d[2] = z - point[2];
  double x0dotr = d[0] * runit[0] + d[1] * runit[1] + d[2] * runit[2];
  c[0] = x0dotr * runit[0];
  c[1] = x0dotr * runit[1];
  c[2] = x0dotr * runit[2];
  a[0] = d[0] - c[0];
  a[1] = d[1] - c[1];
  a[2] = d[2] - c[2];
  b[0] = runit[1] * a[2] - runit[2] * a[1];
  b[1] = runit[2] * a[0] - runit[0] * a[2];
  b[2] = runit[0] * a[1] - runit[1] * a[0];
  disp[0] = a[0] * cosine + b[0] * sine;
  disp[1] = a[1] * cosine + b[1] * sine;
  disp[2] = a[2] * cosine + b[2] * sine;
  x = point[0] + c[0] + disp[0];
  y = point[1] + c[1] + disp[1];
  z = point[2] + c[2] + disp[2];
}
void Region::options(int narg, char **arg)
{
  if (narg < 0) utils::missing_cmd_args(FLERR, "region", error);
  interior = 1;
  scaleflag = 1;
  moveflag = rotateflag = 0;
  openflag = 0;
  for (int i = 0; i < 6; i++) open_faces[i] = 0;
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "units") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "region units", error);
      if (strcmp(arg[iarg + 1], "box") == 0)
        scaleflag = 0;
      else if (strcmp(arg[iarg + 1], "lattice") == 0)
        scaleflag = 1;
      else
        error->all(FLERR, "Illegal region units: {}", arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "side") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "region side", error);
      if (strcmp(arg[iarg + 1], "in") == 0)
        interior = 1;
      else if (strcmp(arg[iarg + 1], "out") == 0)
        interior = 0;
      else
        error->all(FLERR, "Illegal region side: {}", arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "move") == 0) {
      if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, "region move", error);
      if (strcmp(arg[iarg + 1], "NULL") != 0) {
        if (strstr(arg[iarg + 1], "v_") != arg[iarg + 1])
          error->all(FLERR, "Illegal region move x displacement variable: {}", arg[iarg + 1]);
        xstr = utils::strdup(&arg[iarg + 1][2]);
      }
      if (strcmp(arg[iarg + 2], "NULL") != 0) {
        if (strstr(arg[iarg + 2], "v_") != arg[iarg + 2])
          error->all(FLERR, "Illegal region move y displacement variable: {}", arg[iarg + 2]);
        ystr = utils::strdup(&arg[iarg + 2][2]);
      }
      if (strcmp(arg[iarg + 3], "NULL") != 0) {
        if (strstr(arg[iarg + 3], "v_") != arg[iarg + 3])
          error->all(FLERR, "Illegal region move z displacement variable: {}", arg[iarg + 3]);
        zstr = utils::strdup(&arg[iarg + 3][2]);
      }
      moveflag = 1;
      iarg += 4;
    } else if (strcmp(arg[iarg], "rotate") == 0) {
      if (iarg + 8 > narg) utils::missing_cmd_args(FLERR, "region rotate", error);
      if (strstr(arg[iarg + 1], "v_") != arg[iarg + 1]) error->all(FLERR, "Illegal region command");
      tstr = utils::strdup(&arg[iarg + 1][2]);
      point[0] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      point[1] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      point[2] = utils::numeric(FLERR, arg[iarg + 4], false, lmp);
      axis[0] = utils::numeric(FLERR, arg[iarg + 5], false, lmp);
      axis[1] = utils::numeric(FLERR, arg[iarg + 6], false, lmp);
      axis[2] = utils::numeric(FLERR, arg[iarg + 7], false, lmp);
      rotateflag = 1;
      iarg += 8;
    } else if (strcmp(arg[iarg], "open") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "region open", error);
      int iface = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (iface < 1 || iface > 6) error->all(FLERR, "Illegal region open face index: {}", iface);
      open_faces[iface - 1] = 1;
      openflag = 1;
      iarg += 2;
    } else
      error->all(FLERR, "Illegal region command argument: {}", arg[iarg]);
  }
  if ((moveflag || rotateflag) && (strcmp(style, "union") == 0 || strcmp(style, "intersect") == 0))
    error->all(FLERR, "Region union or intersect cannot be dynamic");
  if (scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  } else
    xscale = yscale = zscale = 1.0;
  if (rotateflag) {
    point[0] *= xscale;
    point[1] *= yscale;
    point[2] *= zscale;
  }
  if (rotateflag) {
    double len = sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
    if (len == 0.0) error->all(FLERR, "Region cannot have 0 length rotation vector");
    runit[0] = axis[0] / len;
    runit[1] = axis[1] / len;
    runit[2] = axis[2] / len;
  }
  if (moveflag || rotateflag)
    dynamic = 1;
  else
    dynamic = 0;
}
void Region::point_on_line_segment(double *a, double *b, double *c, double *d)
{
  double ba[3], ca[3];
  MathExtra::sub3(b, a, ba);
  MathExtra::sub3(c, a, ca);
  double t = MathExtra::dot3(ca, ba) / MathExtra::dot3(ba, ba);
  if (t <= 0.0) {
    d[0] = a[0];
    d[1] = a[1];
    d[2] = a[2];
  } else if (t >= 1.0) {
    d[0] = b[0];
    d[1] = b[1];
    d[2] = b[2];
  } else {
    d[0] = a[0] + t * ba[0];
    d[1] = a[1] + t * ba[1];
    d[2] = a[2] + t * ba[2];
  }
}
void Region::set_velocity()
{
  if (vel_timestep == update->ntimestep) return;
  vel_timestep = update->ntimestep;
  if (moveflag) {
    if (update->ntimestep > 0) {
      v[0] = (dx - prev[0]) / update->dt;
      v[1] = (dy - prev[1]) / update->dt;
      v[2] = (dz - prev[2]) / update->dt;
    } else
      v[0] = v[1] = v[2] = 0.0;
    prev[0] = dx;
    prev[1] = dy;
    prev[2] = dz;
  }
  if (rotateflag) {
    rpoint[0] = point[0] + dx;
    rpoint[1] = point[1] + dy;
    rpoint[2] = point[2] + dz;
    if (update->ntimestep > 0) {
      double angvel = (theta - prev[3]) / update->dt;
      omega[0] = angvel * axis[0];
      omega[1] = angvel * axis[1];
      omega[2] = angvel * axis[2];
    } else
      omega[0] = omega[1] = omega[2] = 0.0;
    prev[3] = theta;
  }
  if (varshape) { set_velocity_shape(); }
}
int Region::restart(char *buf, int &n)
{
  int size = *((int *) (&buf[n]));
  n += sizeof(int);
  if ((size <= 0) || (strcmp(&buf[n], id) != 0)) return 0;
  n += size;
  size = *((int *) (&buf[n]));
  n += sizeof(int);
  if ((size <= 0) || (strcmp(&buf[n], style) != 0)) return 0;
  n += size;
  int restart_nreg = *((int *) (&buf[n]));
  n += sizeof(int);
  if (restart_nreg != nregion) return 0;
  memcpy(prev, &buf[n], size_restart * sizeof(double));
  return 1;
}
void Region::reset_vel()
{
  for (int i = 0; i < size_restart; i++) prev[i] = 0;
}
