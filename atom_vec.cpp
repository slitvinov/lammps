#include <set>
#include <map>
#include <unordered_set>
#include "pointers.h"
#include "atom_vec.h"
#include "pointers.h"
#include <vector>
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "fix.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
using namespace LAMMPS_NS;
const std::vector<std::string> AtomVec::default_grow = {
    "id", "type", "mask", "image", "x", "v", "f"};
const std::vector<std::string> AtomVec::default_copy = {"id",    "type", "mask",
                                                        "image", "x",    "v"};
const std::vector<std::string> AtomVec::default_comm = {"x"};
const std::vector<std::string> AtomVec::default_comm_vel = {"x", "v"};
const std::vector<std::string> AtomVec::default_reverse = {"f"};
const std::vector<std::string> AtomVec::default_border = {"id", "type", "mask",
                                                          "x"};
const std::vector<std::string> AtomVec::default_border_vel = {"id", "type",
                                                              "mask", "x", "v"};
const std::vector<std::string> AtomVec::default_exchange = {
    "id", "type", "mask", "image", "x", "v"};
const std::vector<std::string> AtomVec::default_restart = {
    "id", "type", "mask", "image", "x", "v"};
const std::vector<std::string> AtomVec::default_create = {
    "id", "type", "mask", "image", "x", "v"};
const std::vector<std::string> AtomVec::default_data_atom = {};
const std::vector<std::string> AtomVec::default_data_vel = {};
AtomVec::AtomVec(LAMMPS *lmp) : Pointers(lmp) {
  nmax = 0;
  ngrow = 0;
  molecular = Atom::ATOMIC;
  bonds_allow = 0;
  mass_type = dipole_type = PER_ATOM;
  forceclearflag = 0;
  maxexchange = 0;
  bonus_flag = 0;
  size_forward_bonus = size_border_bonus = 0;
  nargcopy = 0;
  argcopy = nullptr;
  tag = nullptr;
  type = mask = nullptr;
  image = nullptr;
  x = v = f = nullptr;
  threads = nullptr;
}
void AtomVec::store_args(int narg, char **arg) {
  nargcopy = narg;
  argcopy = nullptr;
}
void AtomVec::process_args(int narg, char **) {}
void AtomVec::init() {
  deform_vremap = domain->deform_vremap;
  deform_groupbit = domain->deform_groupbit;
  h_rate = domain->h_rate;
}
static constexpr bigint DELTA = 16384;
void AtomVec::grow_nmax() {
  nmax = nmax / DELTA * DELTA;
  nmax += DELTA;
}
static constexpr bigint DELTA_BONUS = 8192;
void AtomVec::grow(int n) {
  int datatype, cols, maxcols;
  void *pdata;
  if (n == 0)
    grow_nmax();
  else
    nmax = MAX(n, nmax);
  atom->nmax = nmax;
  tag = memory->grow(atom->tag, nmax, "atom:tag");
  type = memory->grow(atom->type, nmax, "atom:type");
  mask = memory->grow(atom->mask, nmax, "atom:mask");
  image = memory->grow(atom->image, nmax, "atom:image");
  x = memory->grow(atom->x, nmax, 3, "atom:x");
  v = memory->grow(atom->v, nmax, 3, "atom:v");
  f = memory->grow(atom->f, nmax * comm->nthreads, 3, "atom:f");
  grow_pointers();
}
void AtomVec::copy(int i, int j, int delflag) {
  int m, n, datatype, cols, collength, ncols;
  void *pdata, *plength;
  tag[j] = tag[i];
  type[j] = type[i];
  mask[j] = mask[i];
  image[j] = image[i];
  x[j][0] = x[i][0];
  x[j][1] = x[i][1];
  x[j][2] = x[i][2];
  v[j][0] = v[i][0];
  v[j][1] = v[i][1];
  v[j][2] = v[i][2];
}
void AtomVec::unpack_reverse(int n, int *list, double *buf) {
  int i, j, m, mm, nn, datatype, cols;
  void *pdata;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
  }
}
int AtomVec::pack_border_vel(int n, int *list, double *buf, int pbc_flag,
                             int *pbc) {
  int i, j, m, mm, nn, datatype, cols;
  double dx, dy, dz, dvx, dvy, dvz;
  void *pdata;
  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0] * domain->xprd;
      dy = pbc[1] * domain->yprd;
      dz = pbc[2] * domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
    }
  }
  return m;
}
void AtomVec::unpack_border_vel(int n, int first, double *buf) {
  int i, m, last, mm, nn, datatype, cols;
  void *pdata;
  m = 0;
  last = first + n;
  while (last > nmax)
    grow(0);
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (tagint)ubuf(buf[m++]).i;
    type[i] = (int)ubuf(buf[m++]).i;
    mask[i] = (int)ubuf(buf[m++]).i;
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
  }
}
int AtomVec::size_restart() {
  int i, nn, cols, collength, ncols;
  void *plength;
  int nlocal = atom->nlocal;
  int n = 11 * nlocal;
  if (nrestart) {
    for (nn = 0; nn < nrestart; nn++) {
      cols = mrestart.cols[nn];
      if (cols == 0)
        n += nlocal;
      else if (cols > 0)
        n += cols * nlocal;
      else {
        collength = mrestart.collength[nn];
        plength = mrestart.plength[nn];
        for (i = 0; i < nlocal; i++) {
          if (collength)
            ncols = (*((int ***)plength))[i][collength - 1];
          else
            ncols = (*((int **)plength))[i];
          n += ncols;
        }
      }
    }
  }
  if (bonus_flag)
    n += size_restart_bonus();
  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      for (i = 0; i < nlocal; i++)
        n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);
  return n;
}
void AtomVec::create_atom(int itype, double *coord) {
  int m, n, datatype, cols;
  void *pdata;
  int nlocal = atom->nlocal;
  if (nlocal == nmax)
    grow(0);
  tag[nlocal] = 0;
  type[nlocal] = itype;
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  mask[nlocal] = 1;
  image[nlocal] =
      ((imageint)IMGMAX << IMG2BITS) | ((imageint)IMGMAX << IMGBITS) | IMGMAX;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;
  for (n = 0; n < ncreate; n++) {
    pdata = mcreate.pdata[n];
    datatype = mcreate.datatype[n];
    cols = mcreate.cols[n];
    if (datatype == Atom::DOUBLE) {
      if (cols == 0) {
        double *vec = *((double **)pdata);
        vec[nlocal] = 0.0;
      } else {
        double **array = *((double ***)pdata);
        for (m = 0; m < cols; m++)
          array[nlocal][m] = 0.0;
      }
    } else if (datatype == Atom::INT) {
      if (cols == 0) {
        int *vec = *((int **)pdata);
        vec[nlocal] = 0;
      } else {
        int **array = *((int ***)pdata);
        for (m = 0; m < cols; m++)
          array[nlocal][m] = 0;
      }
    } else if (datatype == Atom::BIGINT) {
      if (cols == 0) {
        bigint *vec = *((bigint **)pdata);
        vec[nlocal] = 0;
      } else {
        bigint **array = *((bigint ***)pdata);
        for (m = 0; m < cols; m++)
          array[nlocal][m] = 0;
      }
    }
  }
  create_atom_post(nlocal);
  atom->nlocal++;
}
void AtomVec::setup_fields() {
  int n, cols;
  ngrow = process_fields(fields_grow, default_grow, &mgrow);
  ncopy = process_fields(fields_copy, default_copy, &mcopy);
  ncomm = process_fields(fields_comm, default_comm, &mcomm);
  ncomm_vel = process_fields(fields_comm_vel, default_comm_vel, &mcomm_vel);
  nreverse = process_fields(fields_reverse, default_reverse, &mreverse);
  nborder = process_fields(fields_border, default_border, &mborder);
  nborder_vel =
      process_fields(fields_border_vel, default_border_vel, &mborder_vel);
  nexchange = process_fields(fields_exchange, default_exchange, &mexchange);
  nrestart = process_fields(fields_restart, default_restart, &mrestart);
  ncreate = process_fields(fields_create, default_create, &mcreate);
  ndata_atom = process_fields(fields_data_atom, default_data_atom, &mdata_atom);
  ndata_vel = process_fields(fields_data_vel, default_data_vel, &mdata_vel);
  init_method(ngrow, &mgrow);
  init_method(ncopy, &mcopy);
  init_method(ncomm, &mcomm);
  init_method(ncomm_vel, &mcomm_vel);
  init_method(nreverse, &mreverse);
  init_method(nborder, &mborder);
  init_method(nborder_vel, &mborder_vel);
  init_method(nexchange, &mexchange);
  init_method(nrestart, &mrestart);
  init_method(ncreate, &mcreate);
  init_method(ndata_atom, &mdata_atom);
  init_method(ndata_vel, &mdata_vel);
  threads = nullptr;
  comm_x_only = 1;
  comm_f_only = 1;
  size_forward = 3;
  size_reverse = 3;
  size_border = 6;
  size_velocity = 3;
  size_data_atom = 0;
  for (n = 0; n < ndata_atom; n++) {
    cols = mdata_atom.cols[n];
    if (atom->peratom[mdata_atom.index[n]].name == "x")
      xcol_data = size_data_atom + 1;
    if (cols == 0)
      size_data_atom++;
    else
      size_data_atom += cols;
  }
  size_data_vel = 0;
  for (n = 0; n < ndata_vel; n++) {
    cols = mdata_vel.cols[n];
    if (cols == 0)
      size_data_vel++;
    else
      size_data_vel += cols;
  }
}
int AtomVec::process_fields(const std::vector<std::string> &words,
                            const std::vector<std::string> &def_words,
                            Method *method) {
  int nfield = words.size();
  int ndef = def_words.size();
  const auto &peratom = atom->peratom;
  const int nperatom = peratom.size();
  method->resize(nfield);
  std::vector<int> &index = method->index;
  int match;
  for (int i = 0; i < nfield; i++) {
    const std::string &field = words[i];
    for (match = 0; match < nperatom; match++)
      if (field == peratom[match].name)
        break;
    index[i] = match;
  }
  return nfield;
}
void AtomVec::init_method(int nfield, Method *method) {
  for (int i = 0; i < nfield; i++) {
    const auto &field = atom->peratom[method->index[i]];
    method->pdata[i] = (void *)field.address;
    method->datatype[i] = field.datatype;
    method->cols[i] = field.cols;
  }
}
void AtomVec::Method::resize(int nfield) {
  pdata.resize(nfield);
  datatype.resize(nfield);
  cols.resize(nfield);
  maxcols.resize(nfield);
  collength.resize(nfield);
  plength.resize(nfield);
  index.resize(nfield);
}
