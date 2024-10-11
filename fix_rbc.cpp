#include "fix_rbc.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"
#include <assert.h>
#include <time.h>
#include <vector>
static const double SMALL = 0.001;
static const double EPSILON = 1.0e-10;
#define DIH_FORCE(i1,i2,i3,i4) \
  do { \
    assert(i1 != i2); \
    assert(i1 != i3); \
    assert(i1 != i4); \
    assert(i2 != i3); \
    assert(i2 != i4); \
    assert(i3 != i4); \
                                                                                            \
    d21x = x[i2][0] - x[i1][0]; \
    d21y = x[i2][1] - x[i1][1]; \
    d21z = x[i2][2] - x[i1][2]; \
    domain->minimum_image(d21x, d21y, d21z); \
                                                                                            \
    d31x = x[i3][0] - x[i1][0]; \
    d31y = x[i3][1] - x[i1][1]; \
    d31z = x[i3][2] - x[i1][2]; \
    domain->minimum_image(d31x, d31y, d31z); \
                                                                                            \
    d32x = x[i3][0] - x[i2][0]; \
    d32y = x[i3][1] - x[i2][1]; \
    d32z = x[i3][2] - x[i2][2]; \
    domain->minimum_image(d32x, d32y, d32z); \
                                                                                            \
    d34x = x[i3][0] - x[i4][0]; \
    d34y = x[i3][1] - x[i4][1]; \
    d34z = x[i3][2] - x[i4][2]; \
    domain->minimum_image(d34x, d34y, d34z); \
                                                                                            \
    d24x = x[i2][0] - x[i4][0]; \
    d24y = x[i2][1] - x[i4][1]; \
    d24z = x[i2][2] - x[i4][2]; \
    domain->minimum_image(d24x, d24y, d24z); \
                                                                                            \
    d14x = x[i1][0] - x[i4][0]; \
    d14y = x[i1][1] - x[i4][1]; \
    d14z = x[i1][2] - x[i4][2]; \
    domain->minimum_image(d14x, d14y, d14z); \
                                                                                            \
    n1x = d21y * d31z - d31y * d21z; \
    n1y = d31x * d21z - d21x * d31z; \
    n1z = d21x * d31y - d31x * d21y; \
    n2x = d34y * d24z - d24y * d34z; \
    n2y = d24x * d34z - d34x * d24z; \
    n2z = d34x * d24y - d24x * d34y; \
    n1 = n1x * n1x + n1y * n1y + n1z * n1z; \
    n2 = n2x * n2x + n2y * n2y + n2z * n2z; \
    nn = sqrt(n1 * n2); \
                                                                                            \
    costheta = (n1x * n2x + n1y * n2y + n1z * n2z) / nn; \
    if (costheta > 1.0) costheta = 1.0; \
    if (costheta < -1.0) costheta = -1.0; \
    sintheta = sqrt(1.0 - costheta * costheta); \
    if (sintheta < SMALL) sintheta = SMALL; \
    mx = -((n1x - n2x) * d14x + (n1y - n2y) * d14y + (n1z - n2z) * d14z); \
    if (mx < 0) sintheta = -sintheta; \
                                                                                            \
    alfa = rbc->kb * rbc->kT * (cos(rbc->theta0) - costheta * sin(rbc->theta0) / sintheta); \
    a11 = -alfa * costheta / n1; \
    a12 = alfa / nn; \
    a22 = -alfa * costheta / n2; \
                                                                                            \
    f1[0] = a11 * (n1y * d32z - n1z * d32y) + a12 * (n2y * d32z - n2z * d32y); \
    f1[1] = a11 * (n1z * d32x - n1x * d32z) + a12 * (n2z * d32x - n2x * d32z); \
    f1[2] = a11 * (n1x * d32y - n1y * d32x) + a12 * (n2x * d32y - n2y * d32x); \
    f2[0] = a11 * (n1z * d31y - n1y * d31z) + a22 * (n2y * d34z - n2z * d34y) + \
        a12 * (n2z * d31y - n2y * d31z + n1y * d34z - n1z * d34y); \
    f2[1] = a11 * (n1x * d31z - n1z * d31x) + a22 * (n2z * d34x - n2x * d34z) + \
        a12 * (n2x * d31z - n2z * d31x + n1z * d34x - n1x * d34z); \
    f2[2] = a11 * (n1y * d31x - n1x * d31y) + a22 * (n2x * d34y - n2y * d34x) + \
        a12 * (n2y * d31x - n2x * d31y + n1x * d34y - n1y * d34x); \
  } while (0)
#define ID2LOCAL(dest,src) \
  do { \
    LAMMPS_NS::tagint o = src + mol * rbc->nv; \
    (dest) = hash_search(tag2i, o, &status); \
    if (status != 0) \
      error->one(FLERR, \
                 "an atom (tag={}, i={}, id={}) has nonlocal connection " \
                 "(tag={}, id={})", \
                 tag, i2, j, o, (src)); \
  } while (0)
struct LAMMPS_NS::FixRBC::RBC {
  int nt, nv;
  int *num_bond;
  int *flip, *halfedge, *next, *tri;
  double a0, area, cq, cutoff, gamc, gamt, kT, ka, kb, kv, l, l0, lambda, q, theta0, volume;
  class NeighList *list;
  int (LAMMPS_NS::Comm::*exchange)(int, double *, double *&);
};
struct LocalAtom {
  int num_bond;
  int *angle;
};
struct Node {
  LAMMPS_NS::tagint key;
  int value;
};
struct Hash {
  size_t M;
  struct Node *nodes;
};
static int hash_ini(size_t, void *, struct Hash *);
static int hash_insert(struct Hash *, LAMMPS_NS::tagint, int);
static int hash_search(struct Hash *, LAMMPS_NS::tagint, int *);
static int local(int i2, LocalAtom *atoms, LAMMPS_NS::tagint tag, double **x, double **v,
                 double **f, LAMMPS_NS::FixRBC::RBC *rbc, LAMMPS_NS::Domain *domain,
                 LAMMPS_NS::Memory *memory, LAMMPS_NS::Error *error, Hash *tag2i)
{
  int h, h0, h1, h3, i0, i1, i3, i4, j, l, nb, status;
  double a11, a12, a22, alfa, costheta, d14x, d14y, d14z, d21x, d21y, d21z, d24x, d24y, d24z, d31x,
      d31y, d31z, d32x, d32y, d32z, d34x, d34y, d34z, delx, dely, delz, fbond, lsq, mx, n1, n1x,
      n1y, n1z, n2, n2x, n2y, n2z, nn, r, rdl, rsq, sintheta, vv, dv[3], f1[3], f2[3];
  LAMMPS_NS::tagint mol;
  j = (tag - 1) % rbc->nv;
  mol = (tag - 1) / rbc->nv;
  nb = atoms->num_bond = rbc->num_bond[j];
  atoms->angle = (int *) memory->smalloc(3 * nb * sizeof(int), "fix/rbc:angle");
  h = h0 = rbc->halfedge[j];
  ID2LOCAL(i0, rbc->tri[rbc->next[rbc->next[rbc->flip[h]]]]);
  ID2LOCAL(i3, rbc->tri[rbc->next[h]]);
  l = 0;
  do {
    h3 = rbc->next[h];
    h1 = rbc->next[h3];
    assert(h == rbc->next[h1]);
    ID2LOCAL(i1, rbc->tri[h1]);
    ID2LOCAL(i4, rbc->tri[rbc->next[rbc->next[rbc->flip[h3]]]]);
    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];
    domain->minimum_image(delx, dely, delz);
    rsq = delx * delx + dely * dely + delz * delz;
    r = sqrt(rsq);
    lsq = rbc->l * rbc->l;
    rdl = sqrt(rsq / lsq);
    fbond = 1.0 / ((1.0 - rdl) * (1.0 - rdl)) + 4.0 * rdl - 1.0;
    fbond *= -0.25 * rbc->kT / rbc->lambda;
    fbond /= r;
    MathExtra::sub3(v[i1], v[i2], dv);
    double del[] = {delx, dely, delz};
    vv = MathExtra::dot3(del, dv) / r;
    fbond -= rbc->gamc * vv / r;
    f[i2][0] -= delx * fbond - rbc->gamt * dv[0];
    f[i2][1] -= dely * fbond - rbc->gamt * dv[1];
    f[i2][2] -= delz * fbond - rbc->gamt * dv[2];
    DIH_FORCE(i1, i2, i3, i0);
    f[i2][0] += f2[0];
    f[i2][1] += f2[1];
    f[i2][2] += f2[2];
    DIH_FORCE(i2, i3, i1, i4);
    f[i2][0] += f1[0];
    f[i2][1] += f1[1];
    f[i2][2] += f1[2];
    atoms->angle[3 * l] = i1;
    atoms->angle[3 * l + 1] = i2;
    atoms->angle[3 * l + 2] = i3;
    i0 = i3;
    i3 = i1;
    h = rbc->flip[h1];
    l++;
  } while (h != h0);
  assert(l == nb);
  return 0;
}
static int area_force(int i2, LocalAtom *atoms, LAMMPS_NS::tagint tag, double **x, double **f,
                      double *center, double *area_volume, LAMMPS_NS::FixRBC::RBC *rbc,
                      LAMMPS_NS::Domain *domain, LAMMPS_NS::Error *error, Hash *index)
{
  int j, i1, i3, k, status;
  double x21[3], x32[3], x31[3], area, c[3], normal[3], coeff_v, coeff_a, totalArea, totalVolume,
      l03, addF[3], tcxmc[3];
  LAMMPS_NS::tagint mol;
  for (j = 0; j < atoms->num_bond; j++) {
    i1 = atoms->angle[3 * j];
    i3 = atoms->angle[3 * j + 2];
    assert(i2 == atoms->angle[3 * j + 1]);
    MathExtra::sub3(x[i2], x[i1], x21);
    MathExtra::sub3(x[i3], x[i2], x32);
    MathExtra::sub3(x[i3], x[i1], x31);
    domain->minimum_image(x21);
    domain->minimum_image(x32);
    domain->minimum_image(x31);
    MathExtra::cross3(x21, x31, normal);
    area = 0.5 * MathExtra::len3(normal);
    mol = (tag - 1) / rbc->nv;
    k = hash_search(index, mol, &status);
    if (status != 0) error->one(FLERR, "hash search failed (mol={})\n", mol);
    totalArea = area_volume[2 * k];
    coeff_a = rbc->q * rbc->cq / (4.0 * pow(area, rbc->q + 2));
    coeff_a +=
        -rbc->ka * rbc->kT * (totalArea - rbc->area) / (4. * area * rbc->area * rbc->l0 * rbc->l0);
    totalVolume = area_volume[2 * k + 1];
    l03 = pow(rbc->l0, 3.);
    coeff_v = -rbc->kv * rbc->kT * (totalVolume - rbc->volume) / (rbc->volume * l03);
    c[0] = (3.0 * x[i2][0] - x21[0] + x32[0]) / 3 - center[3 * k];
    c[1] = (3.0 * x[i2][1] - x21[1] + x32[1]) / 3 - center[3 * k + 1];
    c[2] = (3.0 * x[i2][2] - x21[2] + x32[2]) / 3 - center[3 * k + 2];
    domain->minimum_image(c);
    MathExtra::cross3(normal, x31, addF);
    MathExtra::scale3(-coeff_a, addF);
    f[i2][0] += addF[0];
    f[i2][1] += addF[1];
    f[i2][2] += addF[2];
    MathExtra::cross3(x31, c, tcxmc);
    MathExtra::add3(normal, tcxmc, addF);
    MathExtra::scale3(coeff_v / 6., addF);
    f[i2][0] += addF[0];
    f[i2][1] += addF[1];
    f[i2][2] += addF[2];
  }
  return 0;
}
LAMMPS_NS::FixRBC::FixRBC(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  int a, b, c, i, nvf, status;
  char line[2048], *off_path;
  FILE *off_file;
  size_t nbytes;
  Hash hash;
  void *hash_work;
  arg += 3;
  narg -= 3;
  if (narg == 0) error->all(FLERR, "needs OFF file argument");
  off_path = *arg++;
  narg--;
  rbc = (struct RBC *) memory->smalloc(sizeof *rbc, "fix/rbc:RBC");
  rbc->nt = -1;
  rbc->nv = -1;
  rbc->cutoff = 0.5;
  rbc->a0 = 1000.0;
  rbc->kb = 64.0;
  rbc->kT = 0.0945;
  rbc->l = 1.642599;
  rbc->lambda = 0.00141;
  rbc->theta0 = 6.96;
  rbc->q = 1.0;
  rbc->cq = 1.8;
  rbc->ka = 14193;
  rbc->area = 135;
  rbc->l0 = 0.518;
  rbc->kv = 7352;
  rbc->volume = 94;
  rbc->gamt = 90;
  rbc->gamc = 30;
  rbc->list = NULL;
  rbc->exchange = &LAMMPS_NS::Comm::exchange_variable;
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "rc") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix rbc command: rc");
      rbc->cutoff = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "adpd") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix rbc command: adpd");
      rbc->a0 = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "kb") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix rbc command: kb");
      rbc->kb = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "kT") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix rbc command: kT");
      rbc->kT = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "l") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix rbc command: l");
      rbc->l = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "lambda") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix rbc command: lambda");
      rbc->lambda = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "theta0") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix rbc command: theta0");
      rbc->theta0 = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "q") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix rbc command: q");
      rbc->q = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "cq") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix rbc command: cq");
      rbc->cq = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "ka") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix rbc command: ka");
      rbc->ka = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "area") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix rbc command: area");
      rbc->area = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "l0") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix rbc command: l0");
      rbc->l0 = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "kv") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix rbc command: kv");
      rbc->kv = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "volume") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix rbc command: volume");
      rbc->volume = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "gamt") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix rbc command: gamt");
      rbc->gamt = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "gamc") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix rbc command: gamc");
      rbc->gamc = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "all2all") == 0) {
      rbc->exchange = &LAMMPS_NS::Comm::exchange_variable_all2all;
      if (comm->me == 0) error->warning(FLERR, "fix_rbc: all2all");
      iarg++;
    } else {
      error->all(FLERR, "Unknown fix rbc command {}", arg[iarg]);
    }
  }
  if ((off_file = fopen(off_path, "r")) == NULL) error->all(FLERR, "fail to open '{}'", off_path);
  if (fgets(line, sizeof line, off_file) == NULL) error->all(FLERR, "fail read '{}'", off_path);
  if (strncmp(line, "OFF", 3) != 0) error->all(FLERR, "not an OFF file '{}'", off_path);
  if (fgets(line, sizeof line, off_file) == NULL) error->all(FLERR, "fail read '{}'", off_path);
  if (sscanf(line, "%d %d %*d\n", &rbc->nv, &rbc->nt) != 2)
    error->all(FLERR, "fail read '{}'", off_path);
  rbc->tri = (int *) memory->smalloc(3 * rbc->nt * sizeof(int), "fix/rbc:tri");
  for (i = 0; i < rbc->nv; i++) {
    if (fgets(line, sizeof line, off_file) == NULL) error->all(FLERR, "fail read '{}'", off_path);
  }
  for (i = 0; i < rbc->nt; i++) {
    if (fgets(line, sizeof line, off_file) == NULL) error->all(FLERR, "fail read '{}'", off_path);
    if (sscanf(line, "%d %d %d %d", &nvf, &rbc->tri[3 * i], &rbc->tri[3 * i + 1],
               &rbc->tri[3 * i + 2]) != 4 ||
        nvf != 3)
      error->all(FLERR, "unexpected line '{}' in '{}'", line, off_path);
  }
  for (i = 0; i < 3 * rbc->nt; i++) {
    if (rbc->tri[i] < 0 || rbc->tri[i] >= rbc->nv)
      error->all(FLERR, "invalid vertice refernce {} in '{}', nv = {}", rbc->tri[i], off_path,
                 rbc->nv);
  }
  if (fclose(off_file) != 0) error->all(FLERR, "fail to close '{}'", off_path);
  rbc->num_bond = (int *) memory->smalloc(rbc->nv * sizeof(int), "fix/rbc:num_bond");
  nbytes = 3 * 3 * rbc->nt * sizeof(struct Node);
  hash_work = memory->smalloc(nbytes, "fix/rbc:hash_work");
  hash_ini(nbytes, hash_work, &hash);
  rbc->halfedge = (int *) memory->smalloc(rbc->nv * sizeof(int), "fix/rbc:half_edges/rbc/halfedge");
  rbc->flip = (int *) memory->smalloc(3 * rbc->nt * sizeof(int), "fix/rbc:half_edges/rbc/flip");
  rbc->next = (int *) memory->smalloc(3 * rbc->nt * sizeof(int), "fix/rbc:half_edges/rbc/next");
  for (i = 0; i < rbc->nt; i++) {
    a = rbc->tri[3 * i];
    b = rbc->tri[3 * i + 1];
    c = rbc->tri[3 * i + 2];
    if (hash_insert(&hash, a + b * rbc->nv, 3 * i) != 0) error->one(FLERR, "hash is full");
    if (hash_insert(&hash, b + c * rbc->nv, 3 * i + 1) != 0) error->one(FLERR, "hash is full");
    if (hash_insert(&hash, c + a * rbc->nv, 3 * i + 2) != 0) error->one(FLERR, "hash is full");
    rbc->next[3 * i] = 3 * i + 1;
    rbc->next[3 * i + 1] = 3 * i + 2;
    rbc->next[3 * i + 2] = 3 * i;
  }
  for (i = 0; i < 3 * rbc->nt; i++) rbc->halfedge[rbc->tri[i]] = i;
  for (i = 0; i < rbc->nt; i++) {
    a = rbc->tri[3 * i];
    b = rbc->tri[3 * i + 1];
    c = rbc->tri[3 * i + 2];
    rbc->flip[3 * i] = hash_search(&hash, b + a * rbc->nv, &status);
    rbc->flip[3 * i + 1] = hash_search(&hash, c + b * rbc->nv, &status);
    rbc->flip[3 * i + 2] = hash_search(&hash, a + c * rbc->nv, &status);
  }
  for (i = 0; i < rbc->nv; i++) rbc->num_bond[i] = 0;
  for (i = 0; i < 3 * rbc->nt; i++) rbc->num_bond[rbc->tri[i]]++;
  memory->sfree(hash_work);
}
LAMMPS_NS::FixRBC::~FixRBC()
{
  memory->sfree(rbc->tri);
  memory->sfree(rbc->next);
  memory->sfree(rbc->halfedge);
  memory->sfree(rbc->flip);
  memory->sfree(rbc->num_bond);
  memory->sfree(rbc);
}
int LAMMPS_NS::FixRBC::setmask()
{
  datamask_read = datamask_modify = 0;
  return LAMMPS_NS::FixConst::POST_FORCE;
}
void LAMMPS_NS::FixRBC::post_force(int)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  double **x = atom->x;
  double **f = atom->f;
  double **v = atom->v;
  imageint *image = atom->image;
  std::vector<int> atomlist;
  tagint *tag = atom->tag;
  tagint mol, moli, molj;
  int i, i1, i2, i3, j, jj, jnum, k, natomlist, nin, nmol, nrecv, nsend, status, xbox, ybox, zbox,
      **firstneigh, *jlist, *numneigh;
  double cutsq, delx, dely, delz, fpair, r, rinv, rsq, wd, xtmp, ytmp, ztmp, c[3], normal[3],
      x21[3], x31[3], x32[3], *area_volume, *in, *out;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  void *index_work;
  struct Hash index;
  struct LocalAtom *atoms;
  std::vector<uint64_t> alist;
  std::vector<uint64_t> size;
  std::vector<int> send;
  std::vector<double> center;
  size_t nbytes;
  void *tag2i_work;
  struct Hash tag2i;
  char log_path[FILENAME_MAX];
  static FILE *log_file = NULL;
  struct timespec time;
  natomlist = 0;
  nbytes = 20 * (nall + 1) * sizeof(struct Node);
  tag2i_work = memory->smalloc(nbytes, "fix/rbc:tagi_work");
  hash_ini(nbytes, tag2i_work, &tag2i);
  for (i = nlocal; i < nall; i++) {
    if (mask[i] & groupbit) {
      k = tag[i] - 1;
      if (hash_insert(&tag2i, k, i) != 0) error->one(FLERR, "hash is full");
    }
  }
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      k = tag[i] - 1;
      if (hash_insert(&tag2i, k, i) != 0) error->one(FLERR, "hash is full");
      j = k % rbc->nv;
      natomlist++;
      atomlist.emplace_back(i);
    }
  atoms = (struct LocalAtom *) memory->smalloc(natomlist * sizeof *atoms, "fix/rbc:LocalAtom");
#pragma omp parallel for
  for (int i = 0; i < natomlist; i++) {
    int i2 = atomlist[i];
    local(i2, &atoms[i], tag[i2], x, v, f, rbc, domain, memory, error, &tag2i);
  }
  nbytes = 10 * (natomlist + 1) * sizeof(struct Node);
  index_work = memory->smalloc(nbytes, "fix/rbc:index_work");
  hash_ini(nbytes, index_work, &index);
  nsend = 0;
  nmol = 0;
  for (i = 0; i < natomlist; i++) {
    i2 = atomlist[i];
    mol = (tag[i2] - 1) / rbc->nv;
    j = hash_search(&index, mol, &status);
    if (status != 0 && status != 1) error->one(FLERR, "hash search failed (mol={})\n", mol);
    if (status == 1) {
      j = nmol++;
      if (hash_insert(&index, mol, j) != 0) error->one(FLERR, "hash is full (mol={})\n", mol);
      alist.emplace_back(mol);
      size.emplace_back(0);
      send.emplace_back(1);
      nsend++;
      center.emplace_back(0.0);
      center.emplace_back(0.0);
      center.emplace_back(0.0);
    }
    size[j]++;
    if (size[j] == (uint64_t) rbc->nv) {
      send[j] = 0;
      nsend--;
    }
    xbox = (image[i2] & IMGMASK) - IMGMAX;
    ybox = (image[i2] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (image[i2] >> IMG2BITS) - IMGMAX;
    center[3 * j] += x[i2][0] + xbox * xprd;
    center[3 * j + 1] += x[i2][1] + ybox * yprd;
    center[3 * j + 2] += x[i2][2] + zbox * zprd;
  }
  out = (double *) memory->smalloc(5 * nsend * sizeof *out, "fix/rbc:out");
  for (i = j = 0; i < nmol; i++) {
    if (send[i]) {
      *(uint64_t *) &out[5 * j] = alist[i];
      *(uint64_t *) &out[5 * j + 1] = size[i];
      out[5 * j + 2] = center[3 * i];
      out[5 * j + 3] = center[3 * i + 1];
      out[5 * j + 4] = center[3 * i + 2];
      j++;
    }
  }
  nrecv = (comm->*rbc->exchange)(5 * nsend, out, in);
  assert(nrecv % 5 == 0);
  in += 5 * nsend;
  nin = nrecv / 5 - nsend;
  for (i = 0; i < nin; i++) {
    mol = *(uint64_t *) &in[5 * i];
    if (mol < 0) error->one(FLERR, "mol={} < 0\n", mol);
    j = hash_search(&index, mol, &status);
    if (status != 0 && status != 1) error->one(FLERR, "hash search failed (mol={})\n", mol);
    if (status == 0) {
      if (mol != (tagint) alist[j]) error->one(FLERR, "{} != {}", mol, alist[j]);
      size[j] += *(uint64_t *) &in[5 * i + 1];
      center[3 * j] += in[5 * i + 2];
      center[3 * j + 1] += in[5 * i + 3];
      center[3 * j + 2] += in[5 * i + 4];
    }
  }
  for (i = 0; i < nmol; i++) {
    if (size[i] != (uint64_t) rbc->nv)
      error->one(FLERR, "incomplete cell size[i]={}, rbc->nv={}\n", size[i], rbc->nv);
    center[3 * i] /= size[i];
    center[3 * i + 1] /= size[i];
    center[3 * i + 2] /= size[i];
  }
  area_volume = (double *) memory->smalloc(2 * nmol * sizeof *area_volume, "fix/rbc:area_volume");
  for (i = 0; i < 2 * nmol; i++) area_volume[i] = 0;
  for (i = 0; i < natomlist; i++) {
    i2 = atomlist[i];
    for (j = 0; j < atoms[i].num_bond; j++) {
      i1 = atoms[i].angle[3 * j];
      i3 = atoms[i].angle[3 * j + 2];
      assert(i2 == atoms[i].angle[3 * j + 1]);
      MathExtra::sub3(x[i2], x[i1], x21);
      MathExtra::sub3(x[i3], x[i2], x32);
      MathExtra::sub3(x[i3], x[i1], x31);
      domain->minimum_image(x21);
      domain->minimum_image(x32);
      domain->minimum_image(x31);
      MathExtra::cross3(x21, x31, normal);
      mol = (tag[i2] - 1) / rbc->nv;
      k = hash_search(&index, mol, &status);
      if (status != 0) error->one(FLERR, "hash search failed\n");
      area_volume[2 * k] += (MathExtra::len3(normal) / 2) / 3;
      c[0] = (3.0 * x[i2][0] - x21[0] + x32[0]) / 3 - center[3 * k];
      c[1] = (3.0 * x[i2][1] - x21[1] + x32[1]) / 3 - center[3 * k + 1];
      c[2] = (3.0 * x[i2][2] - x21[2] + x32[2]) / 3 - center[3 * k + 2];
      domain->minimum_image(c);
      area_volume[2 * k + 1] += (MathExtra::dot3(c, normal) / 6) / 3;
    }
  }
  for (i = j = 0; i < nmol; i++)
    if (send[i]) {
      *(uint64_t *) &out[3 * j] = alist[i];
      out[3 * j + 1] = area_volume[2 * i];
      out[3 * j + 2] = area_volume[2 * i + 1];
      j++;
    }
  nrecv = (comm->*rbc->exchange)(3 * nsend, out, in);
  in += 3 * nsend;
  nin = nrecv / 3 - nsend;
  for (i = 0; i < nin; i++) {
    mol = *(uint64_t *) &in[3 * i];
    j = hash_search(&index, mol, &status);
    if (status != 0 && status != 1) error->one(FLERR, "hash search failed (mol={})\n", mol);
    if (status == 0) {
      assert(mol == (tagint) alist[j]);
      area_volume[2 * j] += in[3 * i + 1];
      area_volume[2 * j + 1] += in[3 * i + 2];
    }
  }
#pragma omp parallel for
  for (int i = 0; i < natomlist; i++) {
    int i2 = atomlist[i];
    area_force(i2, &atoms[i], tag[i2], x, f, center.data(), area_volume, rbc, domain, error,
               &index);
  }
  assert(rbc->list);
  if (rbc->a0 > 0) {
    numneigh = rbc->list->numneigh;
    firstneigh = rbc->list->firstneigh;
    cutsq = rbc->cutoff * rbc->cutoff;
    for (k = 0; k < natomlist; k++) {
      i = atomlist[k];
      moli = (tag[i] - 1) / rbc->nv;
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      jlist = firstneigh[i];
      jnum = numneigh[i];
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;
        if (mask[j] & groupbit) {
          molj = (tag[j] - 1) / rbc->nv;
          if (moli != molj) {
            delx = xtmp - x[j][0];
            dely = ytmp - x[j][1];
            delz = ztmp - x[j][2];
            rsq = delx * delx + dely * dely + delz * delz;
            if (rsq < cutsq) {
              const double r = std::sqrt(rsq);
              const double w = 1.0 - r / rbc->cutoff;
              fpair = w * rbc->a0 / (r + 1e-6);
              f[i][0] += delx * fpair;
              f[i][1] += dely * fpair;
              f[i][2] += delz * fpair;
              if (j < nlocal) {
                f[j][0] -= delx * fpair;
                f[j][1] -= dely * fpair;
                f[j][2] -= delz * fpair;
              }
            }
          }
        }
      }
    }
  }
  memory->sfree(tag2i_work);
  for (i = 0; i < natomlist; i++) memory->sfree(atoms[i].angle);
  memory->sfree(atoms);
  memory->sfree(index_work);
  memory->sfree(out);
  memory->sfree(area_volume);
}
void LAMMPS_NS::FixRBC::init_list(int, class NeighList *ptr)
{
  rbc->list = ptr;
}
void LAMMPS_NS::FixRBC::init()
{
  neighbor->add_request(this);
}
static int hash_ini(size_t nbytes, void *memory, struct Hash *hash)
{
  size_t i;
  hash->M = nbytes / sizeof(struct Node);
  hash->nodes = (struct Node *) memory;
  for (i = 0; i < hash->M; i++) hash->nodes[i].key = -1;
  return 0;
}
static int hash_insert(struct Hash *hash, LAMMPS_NS::tagint key, int value)
{
  int x;
  size_t cnt;
  LAMMPS_NS::tagint key0;
  assert(key >= 0);
  x = key % hash->M;
  for (cnt = 0; cnt < hash->M; cnt++) {
    key0 = hash->nodes[x].key;
    if (key0 == -1) {
      hash->nodes[x].key = key;
      hash->nodes[x].value = value;
      return 0;
    } else if (key0 == key) {
      hash->nodes[x].value = value;
      return 0;
    }
    x = (x + 1) % hash->M;
  }
  return -1;
}
static int hash_search(struct Hash *hash, LAMMPS_NS::tagint key, int *status)
{
  int x;
  LAMMPS_NS::tagint key0;
  size_t cnt;
  assert(key >= 0);
  x = key % hash->M;
  for (cnt = 0; cnt < hash->M; cnt++) {
    key0 = hash->nodes[x].key;
    if (key0 == key) {
      *status = 0;
      return hash->nodes[x].value;
    } else if (key0 == -1) {
      *status = 1;
      return -1;
    }
    x = (x + 1) % hash->M;
  }
  *status = 2;
  return -1;
}
