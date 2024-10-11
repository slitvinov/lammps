#include "region_prism.h"
#include "domain.h"
#include "error.h"
#include "math_extra.h"
#include <cstring>
using namespace LAMMPS_NS;
#define BIG 1.0e20
RegPrism::RegPrism(LAMMPS *lmp, int narg, char **arg) : Region(lmp, narg, arg)
{
  options(narg-11,&arg[11]);
  if (strcmp(arg[2],"INF") == 0 || strcmp(arg[2],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[2],"INF") == 0) xlo = -BIG;
    else xlo = domain->boxlo[0];
  } else xlo = xscale*utils::numeric(FLERR,arg[2],false,lmp);
  if (strcmp(arg[3],"INF") == 0 || strcmp(arg[3],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[3],"INF") == 0) xhi = BIG;
    else xhi = domain->boxhi[0];
  } else xhi = xscale*utils::numeric(FLERR,arg[3],false,lmp);
  if (strcmp(arg[4],"INF") == 0 || strcmp(arg[4],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[4],"INF") == 0) ylo = -BIG;
    else ylo = domain->boxlo[1];
  } else ylo = yscale*utils::numeric(FLERR,arg[4],false,lmp);
  if (strcmp(arg[5],"INF") == 0 || strcmp(arg[5],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[5],"INF") == 0) yhi = BIG;
    else yhi = domain->boxhi[1];
  } else yhi = yscale*utils::numeric(FLERR,arg[5],false,lmp);
  if (strcmp(arg[6],"INF") == 0 || strcmp(arg[6],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[6],"INF") == 0) zlo = -BIG;
    else zlo = domain->boxlo[2];
  } else zlo = zscale*utils::numeric(FLERR,arg[6],false,lmp);
  if (strcmp(arg[7],"INF") == 0 || strcmp(arg[7],"EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot use region INF or EDGE when box does not exist");
    if (strcmp(arg[7],"INF") == 0) zhi = BIG;
    else zhi = domain->boxhi[2];
  } else zhi = zscale*utils::numeric(FLERR,arg[7],false,lmp);
  xy = xscale*utils::numeric(FLERR,arg[8],false,lmp);
  xz = xscale*utils::numeric(FLERR,arg[9],false,lmp);
  yz = yscale*utils::numeric(FLERR,arg[10],false,lmp);
  if (xlo >= xhi) error->all(FLERR,"Illegal region prism xlo: {} >= xhi: {}", xlo, xhi);
  if (ylo >= yhi) error->all(FLERR,"Illegal region prism ylo: {} >= yhi: {}", ylo, yhi);
  if (zlo >= zhi) error->all(FLERR,"Illegal region prism zlo: {} >= zhi: {}", zlo ,zhi);
  if (xy != 0.0 && xlo == -BIG && xhi == BIG)
    error->all(FLERR,"Illegal region prism non-zero xy tilt with infinite x size");
  if (xy != 0.0 && ylo == -BIG && yhi == BIG)
    error->all(FLERR,"Illegal region prism non-zero xy tilt with infinite y size");
  if (xz != 0.0 && xlo == -BIG && xhi == BIG)
    error->all(FLERR,"Illegal region prism non-zero xz tilt with infinite x size");
  if (xz != 0.0 && zlo == -BIG && zhi == BIG)
    error->all(FLERR,"Illegal region prism non-zero xz tilt with infinite z size");
  if (yz != 0.0 && ylo == -BIG && yhi == BIG)
    error->all(FLERR,"Illegal region prism non-zero yz tilt with infinite y size");
  if (yz != 0.0 && zlo == -BIG && zhi == BIG)
    error->all(FLERR,"Illegal region prism non-zero yz tilt with infinite z size");
  if (interior) {
    bboxflag = 1;
    extent_xlo = MIN(xlo,xlo+xy);
    extent_xlo = MIN(extent_xlo,extent_xlo+xz);
    extent_ylo = MIN(ylo,ylo+yz);
    extent_zlo = zlo;
    extent_xhi = MAX(xhi,xhi+xy);
    extent_xhi = MAX(extent_xhi,extent_xhi+xz);
    extent_yhi = MAX(yhi,yhi+yz);
    extent_zhi = zhi;
  } else bboxflag = 0;
  cmax = 6;
  contact = new Contact[cmax];
  if (interior) tmax = 3;
  else tmax = 1;
  h[0][0] = xhi - xlo;
  h[0][1] = xy;
  h[0][2] = xz;
  h[1][1] = yhi - ylo;
  h[1][2] = yz;
  h[2][2] = zhi - zlo;
  hinv[0][0] = 1.0/h[0][0];
  hinv[0][1] = -h[0][1] / (h[0][0]*h[1][1]);
  hinv[0][2] = (h[0][1]*h[1][2] - h[0][2]*h[1][1]) / (h[0][0]*h[1][1]*h[2][2]);
  hinv[1][1] = 1.0/h[1][1];
  hinv[1][2] = -h[1][2] / (h[1][1]*h[2][2]);
  hinv[2][2] = 1.0/h[2][2];
  a[0] = xhi-xlo;
  a[1] = 0.0;
  a[2] = 0.0;
  b[0] = xy;
  b[1] = yhi-ylo;
  b[2] = 0.0;
  c[0] = xz;
  c[1] = yz;
  c[2] = zhi-zlo;
  clo[0] = corners[0][0] = xlo;
  clo[1] = corners[0][1] = ylo;
  clo[2] = corners[0][2] = zlo;
  corners[1][0] = xlo + a[0];
  corners[1][1] = ylo + a[1];
  corners[1][2] = zlo + a[2];
  corners[2][0] = xlo + b[0];
  corners[2][1] = ylo + b[1];
  corners[2][2] = zlo + b[2];
  corners[3][0] = xlo + a[0] + b[0];
  corners[3][1] = ylo + a[1] + b[1];
  corners[3][2] = zlo + a[2] + b[2];
  corners[4][0] = xlo + c[0];
  corners[4][1] = ylo + c[1];
  corners[4][2] = zlo + c[2];
  corners[5][0] = xlo + a[0] + c[0];
  corners[5][1] = ylo + a[1] + c[1];
  corners[5][2] = zlo + a[2] + c[2];
  corners[6][0] = xlo + b[0] + c[0];
  corners[6][1] = ylo + b[1] + c[1];
  corners[6][2] = zlo + b[2] + c[2];
  chi[0] = corners[7][0] = xlo + a[0] + b[0] + c[0];
  chi[1] = corners[7][1] = ylo + a[1] + b[1] + c[1];
  chi[2] = corners[7][2] = zlo + a[2] + b[2] + c[2];
  MathExtra::cross3(a,b,face[0]);
  MathExtra::cross3(b,a,face[1]);
  MathExtra::cross3(c,a,face[2]);
  MathExtra::cross3(a,c,face[3]);
  MathExtra::cross3(b,c,face[4]);
  MathExtra::cross3(c,b,face[5]);
  if (openflag) {
    int temp[6];
    for (int i = 0; i < 6; i++)
      temp[i] = open_faces[i];
    open_faces[0] = temp[4];
    open_faces[1] = temp[5];
    open_faces[2] = temp[2];
    open_faces[3] = temp[3];
    open_faces[4] = temp[0];
    open_faces[5] = temp[1];
  }
  for (int i = 0; i < 6; i++) MathExtra::norm3(face[i]);
  tri[0][0] = 0; tri[0][1] = 1; tri[0][2] = 3;
  tri[1][0] = 0; tri[1][1] = 3; tri[1][2] = 2;
  tri[2][0] = 4; tri[2][1] = 7; tri[2][2] = 5;
  tri[3][0] = 4; tri[3][1] = 6; tri[3][2] = 7;
  tri[4][0] = 0; tri[4][1] = 4; tri[4][2] = 5;
  tri[5][0] = 0; tri[5][1] = 5; tri[5][2] = 1;
  tri[6][0] = 2; tri[6][1] = 7; tri[6][2] = 6;
  tri[7][0] = 2; tri[7][1] = 3; tri[7][2] = 7;
  tri[8][0] = 2; tri[8][1] = 6; tri[8][2] = 4;
  tri[9][0] = 2; tri[9][1] = 4; tri[9][2] = 0;
  tri[10][0] = 1; tri[10][1] = 5; tri[10][2] = 7;
  tri[11][0] = 1; tri[11][1] = 7; tri[11][2] = 3;
}
RegPrism::~RegPrism()
{
  delete [] contact;
}
int RegPrism::inside(double x, double y, double z)
{
  double a = hinv[0][0]*(x-xlo) + hinv[0][1]*(y-ylo) + hinv[0][2]*(z-zlo);
  double b = hinv[1][1]*(y-ylo) + hinv[1][2]*(z-zlo);
  double c = hinv[2][2]*(z-zlo);
  if (a >= 0.0 && a <= 1.0 && b >= 0.0 && b <= 1.0 && c >= 0.0 && c <= 1.0)
    return 1;
  return 0;
}
int RegPrism::surface_interior(double *x, double cutoff)
{
  int i;
  double dot;
  double *corner;
  for (i = 0; i < 6; i++) {
    if (i % 2) corner = chi;
    else corner = clo;
    dot = (x[0]-corner[0])*face[i][0] + (x[1]-corner[1])*face[i][1] +
      (x[2]-corner[2])*face[i][2];
    if (dot < 0.0) return 0;
  }
  int n = 0;
  for (i = 0; i < 6; i++) {
    if (open_faces[i]) continue;
    if (i % 2) corner = chi;
    else corner = clo;
    dot = (x[0]-corner[0])*face[i][0] + (x[1]-corner[1])*face[i][1] +
      (x[2]-corner[2])*face[i][2];
    if (dot < cutoff) {
      contact[n].r = dot;
      contact[n].delx = dot*face[i][0];
      contact[n].dely = dot*face[i][1];
      contact[n].delz = dot*face[i][2];
      contact[n].radius = 0;
      contact[n].iwall = i;
      n++;
    }
  }
  return n;
}
int RegPrism::surface_exterior(double *x, double cutoff)
{
  int i;
  double dot;
  double *corner;
  double xp,yp,zp;
  for (i = 0; i < 6; i++) {
    if (i % 2) corner = chi;
    else corner = clo;
    dot = (x[0]-corner[0])*face[i][0] + (x[1]-corner[1])*face[i][1] +
      (x[2]-corner[2])*face[i][2];
    if (dot <= -cutoff) return 0;
  }
  for (i = 0; i < 6; i++) {
    if (i % 2) corner = chi;
    else corner = clo;
    dot = (x[0]-corner[0])*face[i][0] + (x[1]-corner[1])*face[i][1] +
      (x[2]-corner[2])*face[i][2];
    if (dot <= 0.0) break;
  }
  if (i == 6) return 0;
  find_nearest(x,xp,yp,zp);
  add_contact(0,x,xp,yp,zp);
  contact[0].radius = 0;
  contact[0].iwall = 0;
  if (contact[0].r < cutoff) return 1;
  return 0;
}
void RegPrism::find_nearest(double *x, double &xp, double &yp, double &zp)
{
  int i,j,k,iface;
  double xproj[3],xline[3],nearest[3];
  double dot;
  double distsq = BIG;
  for (int itri = 0; itri < 12; itri++) {
    iface = itri/2;
    if (open_faces[iface]) continue;
    i = tri[itri][0];
    j = tri[itri][1];
    k = tri[itri][2];
    dot = (x[0]-corners[i][0])*face[iface][0] +
      (x[1]-corners[i][1])*face[iface][1] +
      (x[2]-corners[i][2])*face[iface][2];
    xproj[0] = x[0] - dot*face[iface][0];
    xproj[1] = x[1] - dot*face[iface][1];
    xproj[2] = x[2] - dot*face[iface][2];
    if (inside_tri(xproj,corners[i],corners[j],corners[k],face[iface])) {
      distsq = closest(x,xproj,nearest,distsq);
    }
    else {
      point_on_line_segment(corners[i],corners[j],xproj,xline);
      distsq = closest(x,xline,nearest,distsq);
      point_on_line_segment(corners[j],corners[k],xproj,xline);
      distsq = closest(x,xline,nearest,distsq);
      point_on_line_segment(corners[i],corners[k],xproj,xline);
      distsq = closest(x,xline,nearest,distsq);
    }
  }
  xp = nearest[0];
  yp = nearest[1];
  zp = nearest[2];
}
int RegPrism::inside_tri(double *x, double *v1, double *v2, double *v3,
                         double *norm)
{
  double edge[3],pvec[3],xproduct[3];
  MathExtra::sub3(v2,v1,edge);
  MathExtra::sub3(x,v1,pvec);
  MathExtra::cross3(edge,pvec,xproduct);
  if (MathExtra::dot3(xproduct,norm) < 0.0) return 0;
  MathExtra::sub3(v3,v2,edge);
  MathExtra::sub3(x,v2,pvec);
  MathExtra::cross3(edge,pvec,xproduct);
  if (MathExtra::dot3(xproduct,norm) < 0.0) return 0;
  MathExtra::sub3(v1,v3,edge);
  MathExtra::sub3(x,v3,pvec);
  MathExtra::cross3(edge,pvec,xproduct);
  if (MathExtra::dot3(xproduct,norm) < 0.0) return 0;
  return 1;
}
double RegPrism::closest(double *x, double *near, double *nearest, double dsq)
{
  double delx = x[0] - near[0];
  double dely = x[1] - near[1];
  double delz = x[2] - near[2];
  double rsq = delx*delx + dely*dely + delz*delz;
  if (rsq >= dsq) return dsq;
  nearest[0] = near[0];
  nearest[1] = near[1];
  nearest[2] = near[2];
  return rsq;
}
