#include "lattice.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "update.h"
#include <cmath>
#include <cstring>
using namespace LAMMPS_NS;
#define BIG 1.0e30
Lattice::Lattice(LAMMPS *lmp, int narg, char **arg) : Pointers(lmp)
{
  nbasis = 0;
  basis = nullptr;
  if (narg < 1) utils::missing_cmd_args(FLERR, "lattice", error);
  if (strcmp(arg[0],"none") == 0) style = NONE;
  else if (strcmp(arg[0],"sc") == 0) style = SC;
  else if (strcmp(arg[0],"bcc") == 0) style = BCC;
  else if (strcmp(arg[0],"fcc") == 0) style = FCC;
  else if (strcmp(arg[0],"hcp") == 0) style = HCP;
  else if (strcmp(arg[0],"diamond") == 0) style = DIAMOND;
  else if (strcmp(arg[0],"sq") == 0) style = SQ;
  else if (strcmp(arg[0],"sq2") == 0) style = SQ2;
  else if (strcmp(arg[0],"hex") == 0) style = HEX;
  else if (strcmp(arg[0],"custom") == 0) style = CUSTOM;
  else error->all(FLERR,"Unknown lattice keyword: {}", arg[0]);
  if (style == NONE) {
    if (narg != 2) error->all(FLERR,"Illegal lattice command: expected 2 arguments but found {}", narg);
    xlattice = ylattice = zlattice = utils::numeric(FLERR,arg[1],false,lmp);
    if (xlattice <= 0.0) error->all(FLERR, "Invalid lattice none argument: {}", arg[1]);
    return;
  }
  int dimension = domain->dimension;
  if (dimension == 2) {
    if (style == SC || style == BCC || style == FCC || style == HCP ||
        style == DIAMOND)
      error->all(FLERR,"Lattice style incompatible with simulation dimension");
  }
  if (dimension == 3) {
    if (style == SQ || style == SQ2 || style == HEX)
      error->all(FLERR,"Lattice style incompatible with simulation dimension");
  }
  if (narg < 2) utils::missing_cmd_args(FLERR, "lattice", error);
  scale = utils::numeric(FLERR,arg[1],false,lmp);
  if (scale <= 0.0) error->all(FLERR, "Invalid lattice {} argument: {}", arg[0], arg[1]);
  if (style == SC) {
    add_basis(0.0,0.0,0.0);
  } else if (style == BCC) {
    add_basis(0.0,0.0,0.0);
    add_basis(0.5,0.5,0.5);
  } else if (style == FCC) {
    add_basis(0.0,0.0,0.0);
    add_basis(0.5,0.5,0.0);
    add_basis(0.5,0.0,0.5);
    add_basis(0.0,0.5,0.5);
  } else if (style == HCP) {
    add_basis(0.0,0.0,0.0);
    add_basis(0.5,0.5,0.0);
    add_basis(0.5,5.0/6.0,0.5);
    add_basis(0.0,1.0/3.0,0.5);
  } else if (style == SQ) {
    add_basis(0.0,0.0,0.0);
  } else if (style == SQ2) {
    add_basis(0.0,0.0,0.0);
    add_basis(0.5,0.5,0.0);
  } else if (style == HEX) {
    add_basis(0.0,0.0,0.0);
    add_basis(0.5,0.5,0.0);
  } else if (style == DIAMOND) {
    add_basis(0.0,0.0,0.0);
    add_basis(0.0,0.5,0.5);
    add_basis(0.5,0.0,0.5);
    add_basis(0.5,0.5,0.0);
    add_basis(0.25,0.25,0.25);
    add_basis(0.25,0.75,0.75);
    add_basis(0.75,0.25,0.75);
    add_basis(0.75,0.75,0.25);
  }
  origin[0] = origin[1] = origin[2] = 0.0;
  orientx[0] = 1; orientx[1] = 0; orientx[2] = 0;
  orienty[0] = 0; orienty[1] = 1; orienty[2] = 0;
  orientz[0] = 0; orientz[1] = 0; orientz[2] = 1;
  int spaceflag = 0;
  a1[0] = 1.0; a1[1] = 0.0; a1[2] = 0.0;
  a2[0] = 0.0; a2[1] = 1.0; a2[2] = 0.0;
  a3[0] = 0.0; a3[1] = 0.0; a3[2] = 1.0;
  if (style == HEX) a2[1] = sqrt(3.0);
  if (style == HCP) {
    a2[1] = sqrt(3.0);
    a3[2] = sqrt(8.0/3.0);
  }
  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"origin") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "lattice origin", error);
      origin[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      origin[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      origin[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      if (origin[0] < 0.0 || origin[0] >= 1.0)
        error->all(FLERR, "Invalid lattice origin argument: {}", origin[0]);
      if (origin[1] < 0.0 || origin[1] >= 1.0)
        error->all(FLERR, "Invalid lattice origin argument: {}", origin[1]);
      if (origin[2] < 0.0 || origin[2] >= 1.0)
        error->all(FLERR, "Invalid lattice origin argument: {}", origin[2]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"orient") == 0) {
      if (iarg+5 > narg) utils::missing_cmd_args(FLERR, "lattice orient", error);
      int dim = -1;
      if (strcmp(arg[iarg+1],"x") == 0) dim = 0;
      else if (strcmp(arg[iarg+1],"y") == 0) dim = 1;
      else if (strcmp(arg[iarg+1],"z") == 0) dim = 2;
      else error->all(FLERR,"Unknown lattice orient argument: {}", arg[iarg+1]);
      int *orient = nullptr;
      if (dim == 0) orient = orientx;
      else if (dim == 1) orient = orienty;
      else if (dim == 2) orient = orientz;
      orient[0] = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
      orient[1] = utils::inumeric(FLERR,arg[iarg+3],false,lmp);
      orient[2] = utils::inumeric(FLERR,arg[iarg+4],false,lmp);
      iarg += 5;
    } else if (strcmp(arg[iarg],"spacing") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "lattice spacing", error);
      spaceflag = 1;
      xlattice = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      ylattice = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      zlattice = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      iarg += 4;
    } else if (strcmp(arg[iarg],"a1") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "lattice a1", error);
      if (style != CUSTOM)
        error->all(FLERR,
                   "Invalid a1 option in lattice command for non-custom style");
      a1[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      a1[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      a1[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      iarg += 4;
    } else if (strcmp(arg[iarg],"a2") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "lattice a2", error);
      if (style != CUSTOM)
        error->all(FLERR,
                   "Invalid a2 option in lattice command for non-custom style");
      a2[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      a2[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      a2[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      iarg += 4;
    } else if (strcmp(arg[iarg],"a3") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "lattice a3", error);
      if (style != CUSTOM)
        error->all(FLERR,
                   "Invalid a3 option in lattice command for non-custom style");
      a3[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      a3[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      a3[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      iarg += 4;
    } else if (strcmp(arg[iarg],"basis") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "lattice basis", error);
      if (style != CUSTOM)
        error->all(FLERR,
                   "Invalid basis option in lattice command for non-custom style");
      double x = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      double y = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      double z = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      if (x < 0.0 || x >= 1.0)
        error->all(FLERR, "Invalid lattice basis argument: {}", x);
      if (y < 0.0 || y >= 1.0)
        error->all(FLERR, "Invalid lattice basis argument: {}", y);
      if (z < 0.0 || z >= 1.0)
        error->all(FLERR, "Invalid lattice basis argument: {}", z);
      add_basis(x,y,z);
      iarg += 4;
    } else error->all(FLERR,"Unknown lattice keyword: {}", arg[iarg]);
  }
  if (nbasis == 0) error->all(FLERR,"No basis atoms in lattice");
  if (!orthogonal())
    error->all(FLERR,"Lattice orient vectors are not orthogonal");
  if (!right_handed())
    error->all(FLERR,"Lattice orient vectors are not right-handed");
  if (collinear())
    error->all(FLERR,"Lattice primitive vectors are collinear");
  if (dimension == 2) {
    if (origin[2] != 0.0)
      error->all(FLERR,
                 "Lattice settings are not compatible with 2d simulation");
    if (orientx[2] != 0 || orienty[2] != 0 ||
        orientz[0] != 0 || orientz[1] != 0)
      error->all(FLERR,
                 "Lattice settings are not compatible with 2d simulation");
    if (a1[2] != 0.0 || a2[2] != 0.0 || a3[0] != 0.0 || a3[1] != 0.0)
      error->all(FLERR,
                 "Lattice settings are not compatible with 2d simulation");
  }
  if (spaceflag) {
    if (xlattice <= 0.0 || ylattice <= 0.0 || zlattice <= 0.0)
      error->all(FLERR,"Lattice spacings are invalid");
  }
  if (strcmp(update->unit_style,"lj") == 0) {
    double vec[3];
    cross(a2,a3,vec);
    double volume = dot(a1,vec);
    scale = pow(nbasis/volume/scale,1.0/dimension);
  }
  setup_transform();
  if (spaceflag == 0) {
    double xmin,ymin,zmin,xmax,ymax,zmax;
    xmin = ymin = zmin = BIG;
    xmax = ymax = zmax = -BIG;
    xlattice = ylattice = zlattice = 0.0;
    bbox(0,0.0,0.0,0.0,xmin,ymin,zmin,xmax,ymax,zmax);
    bbox(0,1.0,0.0,0.0,xmin,ymin,zmin,xmax,ymax,zmax);
    bbox(0,0.0,1.0,0.0,xmin,ymin,zmin,xmax,ymax,zmax);
    bbox(0,1.0,1.0,0.0,xmin,ymin,zmin,xmax,ymax,zmax);
    bbox(0,0.0,0.0,1.0,xmin,ymin,zmin,xmax,ymax,zmax);
    bbox(0,1.0,0.0,1.0,xmin,ymin,zmin,xmax,ymax,zmax);
    bbox(0,0.0,1.0,1.0,xmin,ymin,zmin,xmax,ymax,zmax);
    bbox(0,1.0,1.0,1.0,xmin,ymin,zmin,xmax,ymax,zmax);
    xlattice = xmax - xmin;
    ylattice = ymax - ymin;
    zlattice = zmax - zmin;
  } else {
    xlattice *= scale;
    ylattice *= scale;
    zlattice *= scale;
  }
  if (comm->me == 0)
    utils::logmesg(lmp,"Lattice spacing in x,y,z = {:.8} {:.8} {:.8}\n",
                   xlattice,ylattice,zlattice);
}
Lattice::~Lattice()
{
  memory->destroy(basis);
}
int Lattice::orthogonal()
{
  if (orientx[0]*orienty[0] + orientx[1]*orienty[1] +
      orientx[2]*orienty[2]) return 0;
  if (orienty[0]*orientz[0] + orienty[1]*orientz[1] +
      orienty[2]*orientz[2]) return 0;
  if (orientx[0]*orientz[0] + orientx[1]*orientz[1] +
      orientx[2]*orientz[2]) return 0;
  return 1;
}
int Lattice::right_handed()
{
  int xy0 = orientx[1]*orienty[2] - orientx[2]*orienty[1];
  int xy1 = orientx[2]*orienty[0] - orientx[0]*orienty[2];
  int xy2 = orientx[0]*orienty[1] - orientx[1]*orienty[0];
  if (xy0*orientz[0] + xy1*orientz[1] + xy2*orientz[2] <= 0) return 0;
  return 1;
}
int Lattice::collinear()
{
  double vec[3];
  cross(a1,a2,vec);
  if (dot(vec,vec) == 0.0) return 1;
  cross(a2,a3,vec);
  if (dot(vec,vec) == 0.0) return 1;
  cross(a1,a3,vec);
  if (dot(vec,vec) == 0.0) return 1;
  return 0;
}
void Lattice::setup_transform()
{
  double length;
  primitive[0][0] = a1[0];
  primitive[1][0] = a1[1];
  primitive[2][0] = a1[2];
  primitive[0][1] = a2[0];
  primitive[1][1] = a2[1];
  primitive[2][1] = a2[2];
  primitive[0][2] = a3[0];
  primitive[1][2] = a3[1];
  primitive[2][2] = a3[2];
  double determinant = primitive[0][0]*primitive[1][1]*primitive[2][2] +
    primitive[0][1]*primitive[1][2]*primitive[2][0] +
    primitive[0][2]*primitive[1][0]*primitive[2][1] -
    primitive[0][0]*primitive[1][2]*primitive[2][1] -
    primitive[0][1]*primitive[1][0]*primitive[2][2] -
    primitive[0][2]*primitive[1][1]*primitive[2][0];
  if (determinant == 0.0)
    error->all(FLERR,"Degenerate lattice primitive vectors");
  priminv[0][0] = (primitive[1][1]*primitive[2][2] -
                   primitive[1][2]*primitive[2][1]) / determinant;
  priminv[1][0] = (primitive[1][2]*primitive[2][0] -
                   primitive[1][0]*primitive[2][2]) / determinant;
  priminv[2][0] = (primitive[1][0]*primitive[2][1] -
                   primitive[1][1]*primitive[2][0]) / determinant;
  priminv[0][1] = (primitive[0][2]*primitive[2][1] -
                   primitive[0][1]*primitive[2][2]) / determinant;
  priminv[1][1] = (primitive[0][0]*primitive[2][2] -
                   primitive[0][2]*primitive[2][0]) / determinant;
  priminv[2][1] = (primitive[0][1]*primitive[2][0] -
                   primitive[0][0]*primitive[2][1]) / determinant;
  priminv[0][2] = (primitive[0][1]*primitive[1][2] -
                   primitive[0][2]*primitive[1][1]) / determinant;
  priminv[1][2] = (primitive[0][2]*primitive[1][0] -
                   primitive[0][0]*primitive[1][2]) / determinant;
  priminv[2][2] = (primitive[0][0]*primitive[1][1] -
                   primitive[0][1]*primitive[1][0]) / determinant;
  int lensq = orientx[0]*orientx[0] + orientx[1]*orientx[1] +
    orientx[2]*orientx[2];
  length = sqrt((double) lensq);
  if (length == 0.0) error->all(FLERR,"Zero-length lattice orient vector");
  rotaterow[0][0] = orientx[0] / length;
  rotaterow[0][1] = orientx[1] / length;
  rotaterow[0][2] = orientx[2] / length;
  lensq = orienty[0]*orienty[0] + orienty[1]*orienty[1] +
    orienty[2]*orienty[2];
  length = sqrt((double) lensq);
  if (length == 0.0) error->all(FLERR,"Zero-length lattice orient vector");
  rotaterow[1][0] = orienty[0] / length;
  rotaterow[1][1] = orienty[1] / length;
  rotaterow[1][2] = orienty[2] / length;
  lensq = orientz[0]*orientz[0] + orientz[1]*orientz[1] +
    orientz[2]*orientz[2];
  length = sqrt((double) lensq);
  if (length == 0.0) error->all(FLERR,"Zero-length lattice orient vector");
  rotaterow[2][0] = orientz[0] / length;
  rotaterow[2][1] = orientz[1] / length;
  rotaterow[2][2] = orientz[2] / length;
  rotatecol[0][0] = rotaterow[0][0];
  rotatecol[1][0] = rotaterow[0][1];
  rotatecol[2][0] = rotaterow[0][2];
  rotatecol[0][1] = rotaterow[1][0];
  rotatecol[1][1] = rotaterow[1][1];
  rotatecol[2][1] = rotaterow[1][2];
  rotatecol[0][2] = rotaterow[2][0];
  rotatecol[1][2] = rotaterow[2][1];
  rotatecol[2][2] = rotaterow[2][2];
}
void Lattice::lattice2box(double &x, double &y, double &z)
{
  double x1 = primitive[0][0]*x + primitive[0][1]*y + primitive[0][2]*z;
  double y1 = primitive[1][0]*x + primitive[1][1]*y + primitive[1][2]*z;
  double z1 = primitive[2][0]*x + primitive[2][1]*y + primitive[2][2]*z;
  x1 *= scale;
  y1 *= scale;
  z1 *= scale;
  double xnew = rotaterow[0][0]*x1 + rotaterow[0][1]*y1 + rotaterow[0][2]*z1;
  double ynew = rotaterow[1][0]*x1 + rotaterow[1][1]*y1 + rotaterow[1][2]*z1;
  double znew = rotaterow[2][0]*x1 + rotaterow[2][1]*y1 + rotaterow[2][2]*z1;
  x = xnew + xlattice*origin[0];
  y = ynew + ylattice*origin[1];
  z = znew + zlattice*origin[2];
}
void Lattice::box2lattice(double &x, double &y, double &z)
{
  x -= xlattice*origin[0];
  y -= ylattice*origin[1];
  z -= zlattice*origin[2];
  double x1 = rotatecol[0][0]*x + rotatecol[0][1]*y + rotatecol[0][2]*z;
  double y1 = rotatecol[1][0]*x + rotatecol[1][1]*y + rotatecol[1][2]*z;
  double z1 = rotatecol[2][0]*x + rotatecol[2][1]*y + rotatecol[2][2]*z;
  x1 /= scale;
  y1 /= scale;
  z1 /= scale;
  x = priminv[0][0]*x1 + priminv[0][1]*y1 + priminv[0][2]*z1;
  y = priminv[1][0]*x1 + priminv[1][1]*y1 + priminv[1][2]*z1;
  z = priminv[2][0]*x1 + priminv[2][1]*y1 + priminv[2][2]*z1;
}
void Lattice::add_basis(double x, double y, double z)
{
  memory->grow(basis,nbasis+1,3,"lattice:basis");
  basis[nbasis][0] = x;
  basis[nbasis][1] = y;
  basis[nbasis][2] = z;
  nbasis++;
}
double Lattice::dot(double *x, double *y)
{
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}
void Lattice::cross(double *x, double *y, double *z)
{
  z[0] = x[1]*y[2] - x[2]*y[1];
  z[1] = x[2]*y[0] - x[0]*y[2];
  z[2] = x[0]*y[1] - x[1]*y[0];
}
void Lattice::bbox(int flag, double x, double y, double z,
                   double &xmin, double &ymin, double &zmin,
                   double &xmax, double &ymax, double &zmax)
{
  if (flag == 0) lattice2box(x,y,z);
  else box2lattice(x,y,z);
  xmin = MIN(x,xmin); ymin = MIN(y,ymin); zmin = MIN(z,zmin);
  xmax = MAX(x,xmax); ymax = MAX(y,ymax); zmax = MAX(z,zmax);
}
