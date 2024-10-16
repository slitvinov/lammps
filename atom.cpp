#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "library.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "style_atom.h"
#include "update.h"
#include <algorithm>
#include <cstring>
using namespace LAMMPS_NS;
using namespace MathConst;
#define DELTA 1
#define EPSILON 1.0e-6
#define MAXLINE 256
template <typename T> static AtomVec *avec_creator(LAMMPS *_lmp) {
  return new T(_lmp);
}
Atom::Atom(LAMMPS *_lmp) : Pointers(_lmp) {
  natoms = 0;
  nlocal = nghost = nmax = 0;
  ntypes = 0;
  firstgroupname = nullptr;
  sortfreq = 1000;
  nextsort = 0;
  userbinsize = 0.0;
  maxbin = maxnext = 0;
  binhead = nullptr;
  next = permute = nullptr;
  tag = nullptr;
  type = mask = nullptr;
  image = nullptr;
  x = v = f = nullptr;
  q = nullptr;
  mu = nullptr;
  omega = angmom = torque = nullptr;
  radius = rmass = nullptr;
  quat = nullptr;
  temperature = nullptr;
  heatflow = nullptr;
  vfrac = s0 = nullptr;
  x0 = nullptr;
  sp = fm = fm_long = nullptr;
  spin = nullptr;
  eradius = ervel = erforce = nullptr;
  ervelforce = nullptr;
  cs = csforce = vforce = nullptr;
  etag = nullptr;
  id5p = nullptr;
  cc = cc_flux = nullptr;
  edpd_temp = edpd_flux = edpd_cv = nullptr;
  contact_radius = nullptr;
  smd_data_9 = nullptr;
  smd_stress = nullptr;
  eff_plastic_strain = nullptr;
  eff_plastic_strain_rate = nullptr;
  damage = nullptr;
  rho = drho = esph = desph = cv = nullptr;
  vest = nullptr;
  types_style = NUMERIC;
  nivector = ndvector = niarray = ndarray = 0;
  ivector = nullptr;
  dvector = nullptr;
  iarray = nullptr;
  darray = nullptr;
  icols = dcols = nullptr;
  ivname = dvname = ianame = daname = nullptr;
  set_atomflag_defaults();
  peratom_create();
  mass = nullptr;
  mass_setflag = nullptr;
  nextra_grow = nextra_restart = nextra_border = 0;
  extra_grow = extra_restart = extra_border = nullptr;
  nextra_grow_max = nextra_restart_max = nextra_border_max = 0;
  nextra_store = 0;
  extra = nullptr;
  tag_enable = 1;
  map_style = map_user = MAP_NONE;
  map_tag_max = -1;
  map_maxarray = map_nhash = map_nbucket = -1;
  max_same = 0;
  sametag = nullptr;
  map_array = nullptr;
  map_bucket = nullptr;
  map_hash = nullptr;
  unique_tags = nullptr;
  reset_image_flag[0] = reset_image_flag[1] = reset_image_flag[2] = false;
  atom_style = nullptr;
  avec = nullptr;
  avec_map = new AtomVecCreatorMap();
#define ATOM_CLASS
#define AtomStyle(key, Class) (*avec_map)[#key] = &avec_creator<Class>;
#include "style_atom.h"
#undef AtomStyle
#undef ATOM_CLASS
}
Atom::~Atom() {
  delete[] atom_style;
  delete avec;
  delete avec_map;
  delete[] firstgroupname;
  memory->destroy(binhead);
  memory->destroy(next);
  memory->destroy(permute);
  memory->destroy(tag);
  memory->destroy(type);
  memory->destroy(mask);
  memory->destroy(image);
  memory->destroy(x);
  memory->destroy(v);
  memory->destroy(f);
  memory->sfree(ivname);
  memory->sfree(dvname);
  memory->sfree(ianame);
  memory->sfree(daname);
  memory->sfree(ivector);
  memory->sfree(dvector);
  memory->sfree(iarray);
  memory->sfree(darray);
  memory->sfree(icols);
  memory->sfree(dcols);
  delete[] mass;
  delete[] mass_setflag;
  memory->destroy(extra_grow);
  memory->destroy(extra_restart);
  memory->destroy(extra_border);
  memory->destroy(extra);
  Atom::map_delete();
  delete unique_tags;
}
void Atom::peratom_create() {
  peratom.clear();
  int tagintsize = INT;
  if (sizeof(tagint) == 8)
    tagintsize = BIGINT;
  int imageintsize = INT;
  if (sizeof(imageint) == 8)
    imageintsize = BIGINT;
  add_peratom("id", &tag, tagintsize, 0);
  add_peratom("type", &type, INT, 0);
  add_peratom("mask", &mask, INT, 0);
  add_peratom("image", &image, imageintsize, 0);
  add_peratom("x", &x, DOUBLE, 3);
  add_peratom("v", &v, DOUBLE, 3);
  add_peratom("f", &f, DOUBLE, 3, 1);
  add_peratom("rmass", &rmass, DOUBLE, 0);
  add_peratom("q", &q, DOUBLE, 0);
  add_peratom("mu", &mu, DOUBLE, 4);
  add_peratom("mu3", &mu, DOUBLE, 3);
  add_peratom("radius", &radius, DOUBLE, 0);
  add_peratom("omega", &omega, DOUBLE, 3);
  add_peratom("torque", &torque, DOUBLE, 3, 1);
  add_peratom("angmom", &angmom, DOUBLE, 3);
  add_peratom("temperature", &temperature, DOUBLE, 0);
  add_peratom("heatflow", &heatflow, DOUBLE, 0);
  add_peratom("quat", &quat, DOUBLE, 4);
  add_peratom("id5p", &id5p, tagintsize, 0);
  add_peratom("edpd_cv", &edpd_cv, DOUBLE, 0);
  add_peratom("edpd_temp", &edpd_temp, DOUBLE, 0);
  add_peratom("vest_temp", &vest_temp, DOUBLE, 0);
  add_peratom("edpd_flux", &edpd_flux, DOUBLE, 0, 1);
  add_peratom("cc", &cc, DOUBLE, 1);
  add_peratom("cc_flux", &cc_flux, DOUBLE, 1, 1);
  add_peratom("rho", &rho, DOUBLE, 0);
  add_peratom("drho", &drho, DOUBLE, 0, 1);
  add_peratom("esph", &esph, DOUBLE, 0);
  add_peratom("desph", &desph, DOUBLE, 0, 1);
  add_peratom("vest", &vest, DOUBLE, 3);
  add_peratom("cv", &cv, DOUBLE, 0);
  add_peratom("contact_radius", &contact_radius, DOUBLE, 0);
  add_peratom("smd_data_9", &smd_data_9, DOUBLE, 1);
  add_peratom("smd_stress", &smd_stress, DOUBLE, 1);
  add_peratom("eff_plastic_strain", &eff_plastic_strain, DOUBLE, 0);
  add_peratom("eff_plastic_strain_rate", &eff_plastic_strain_rate, DOUBLE, 0);
  add_peratom("damage", &damage, DOUBLE, 0);
}
void Atom::add_peratom(const std::string &name, void *address, int datatype,
                       int cols, int threadflag) {
  PerAtom item = {name,     address, nullptr, nullptr,
                  datatype, cols,    0,       threadflag};
  peratom.push_back(item);
}
void Atom::set_atomflag_defaults() {
  quat_flag = 0;
  peri_flag = electron_flag = 0;
  wavepacket_flag = sph_flag = 0;
  q_flag = mu_flag = 0;
  rmass_flag = radius_flag = omega_flag = torque_flag = angmom_flag = 0;
  temperature_flag = heatflow_flag = 0;
  vfrac_flag = spin_flag = eradius_flag = ervel_flag = erforce_flag = 0;
  cs_flag = csforce_flag = vforce_flag = ervelforce_flag = etag_flag = 0;
  rho_flag = esph_flag = cv_flag = vest_flag = 0;
  dpd_flag = edpd_flag = tdpd_flag = 0;
}
void Atom::create_avec(const std::string &style, int narg, char **arg,
                       int trysuffix) {
  delete[] atom_style;
  if (avec)
    delete avec;
  atom_style = nullptr;
  avec = nullptr;
  set_atomflag_defaults();
  avec = new_avec(style);
  avec->store_args(narg, arg);
  avec->process_args(narg, arg);
  avec->grow(1);
  atom_style = utils::strdup(style);
}
AtomVec *Atom::new_avec(const std::string &style) {
  if (avec_map->find(style) != avec_map->end()) {
    AtomVecCreator &avec_creator = (*avec_map)[style];
    return avec_creator(lmp);
  }
  return nullptr;
}
void Atom::init() {
  if (nextra_store) {
    memory->destroy(extra);
    extra = nullptr;
    nextra_store = 0;
  }
  check_mass(FLERR);
  if (firstgroupname) {
    firstgroup = group->find(firstgroupname);
    if (firstgroup < 0)
      error->all(FLERR, "Could not find atom_modify first group ID {}",
                 firstgroupname);
  } else
    firstgroup = -1;
  avec->init();
}
void Atom::setup() {
  if (sortfreq > 0)
    setup_sort_bins();
}
void Atom::tag_check() {
  tagint min = MAXTAGINT;
  tagint max = 0;
  for (int i = 0; i < nlocal; i++) {
    min = MIN(min, tag[i]);
    max = MAX(max, tag[i]);
  }
  tagint minall, maxall;
  MPI_Allreduce(&min, &minall, 1, MPI_LMP_TAGINT, MPI_MIN, world);
  MPI_Allreduce(&max, &maxall, 1, MPI_LMP_TAGINT, MPI_MAX, world);
  if (minall < 0)
    error->all(FLERR, "One or more Atom IDs are negative");
  if (maxall >= MAXTAGINT)
    error->all(FLERR, "One or more atom IDs are too big");
  if (maxall > 0 && minall == 0)
    error->all(FLERR, "One or more atom IDs are zero");
  if (maxall > 0 && tag_enable == 0)
    error->all(FLERR, "Non-zero atom IDs with atom_modify id = no");
  if (maxall == 0 && natoms && tag_enable)
    error->all(FLERR, "All atom IDs = 0 but atom_modify id = yes");
  if (tag_enable && maxall < natoms)
    error->all(FLERR, "Duplicate atom IDs exist");
}
void Atom::tag_extend() {
  tagint maxtag = 0;
  for (int i = 0; i < nlocal; i++)
    maxtag = MAX(maxtag, tag[i]);
  tagint maxtag_all;
  MPI_Allreduce(&maxtag, &maxtag_all, 1, MPI_LMP_TAGINT, MPI_MAX, world);
  bigint notag = 0;
  for (int i = 0; i < nlocal; i++)
    if (tag[i] == 0)
      notag++;
  bigint notag_total;
  MPI_Allreduce(&notag, &notag_total, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if (notag_total >= MAXTAGINT)
    error->all(FLERR, "New atom IDs exceed maximum allowed ID {}", MAXTAGINT);
  bigint notag_sum;
  MPI_Scan(&notag, &notag_sum, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  tagint itag = maxtag_all + notag_sum - notag + 1;
  for (int i = 0; i < nlocal; i++)
    if (tag[i] == 0)
      tag[i] = itag++;
}
int Atom::tag_consecutive() {
  tagint idmin = MAXTAGINT;
  tagint idmax = 0;
  for (int i = 0; i < nlocal; i++) {
    idmin = MIN(idmin, tag[i]);
    idmax = MAX(idmax, tag[i]);
  }
  tagint idminall, idmaxall;
  MPI_Allreduce(&idmin, &idminall, 1, MPI_LMP_TAGINT, MPI_MIN, world);
  MPI_Allreduce(&idmax, &idmaxall, 1, MPI_LMP_TAGINT, MPI_MAX, world);
  if (idminall != 1 || idmaxall != natoms)
    return 0;
  return 1;
}
void Atom::allocate_type_arrays() {
  if (avec->mass_type == AtomVec::PER_TYPE) {
    mass = new double[ntypes + 1];
    mass_setflag = new int[ntypes + 1];
    for (int itype = 1; itype <= ntypes; itype++)
      mass_setflag[itype] = 0;
  }
}
void Atom::set_mass(const char *file, int line, int, char **arg) {
  if (mass == nullptr)
    error->all(file, line, "Cannot set per-type atom mass for atom style {}",
               atom_style);
  const std::string str = arg[0];
  int lo, hi;
  utils::bounds(file, line, str, 1, ntypes, lo, hi, error);
  if ((lo < 1) || (hi > ntypes))
    error->all(file, line, "Invalid atom type {} for atom mass", str);
  const double value = utils::numeric(FLERR, arg[1], false, lmp);
  if (value <= 0.0)
    error->all(file, line, "Invalid atom mass value {} for type {}", value,
               str);
  for (int itype = lo; itype <= hi; itype++) {
    mass[itype] = value;
    mass_setflag[itype] = 1;
  }
}
void Atom::check_mass(const char *file, int line) {
  if (mass == nullptr)
    return;
  if (rmass_flag)
    return;
  for (int itype = 1; itype <= ntypes; itype++)
    if (mass_setflag[itype] == 0)
      error->all(file, line,
                 "Not all per-type masses are set. Type {} is missing.", itype);
}
void Atom::first_reorder() {
  if (nlocal == nmax)
    avec->grow(0);
  int bitmask = group->bitmask[firstgroup];
  nfirst = 0;
  while (nfirst < nlocal && mask[nfirst] & bitmask)
    nfirst++;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & bitmask && i > nfirst) {
      avec->copy(i, nlocal, 0);
      avec->copy(nfirst, i, 0);
      avec->copy(nlocal, nfirst, 0);
      while (nfirst < nlocal && mask[nfirst] & bitmask)
        nfirst++;
    }
  }
}
void Atom::sort() {
  int i, m, n, ix, iy, iz, ibin, empty;
  nextsort = (update->ntimestep / sortfreq) * sortfreq + sortfreq;
  if (domain->box_change)
    setup_sort_bins();
  if (nbins == 1)
    return;
  if (nlocal > maxnext) {
    memory->destroy(next);
    memory->destroy(permute);
    maxnext = atom->nmax;
    memory->create(next, maxnext, "atom:next");
    memory->create(permute, maxnext, "atom:permute");
  }
  if (nlocal == nmax)
    avec->grow(0);
  for (i = 0; i < nbins; i++)
    binhead[i] = -1;
  for (i = nlocal - 1; i >= 0; i--) {
    ix = static_cast<int>((x[i][0] - bboxlo[0]) * bininvx);
    iy = static_cast<int>((x[i][1] - bboxlo[1]) * bininvy);
    iz = static_cast<int>((x[i][2] - bboxlo[2]) * bininvz);
    ix = MAX(ix, 0);
    iy = MAX(iy, 0);
    iz = MAX(iz, 0);
    ix = MIN(ix, nbinx - 1);
    iy = MIN(iy, nbiny - 1);
    iz = MIN(iz, nbinz - 1);
    ibin = iz * nbiny * nbinx + iy * nbinx + ix;
    next[i] = binhead[ibin];
    binhead[ibin] = i;
  }
  n = 0;
  for (m = 0; m < nbins; m++) {
    i = binhead[m];
    while (i >= 0) {
      permute[n++] = i;
      i = next[i];
    }
  }
  int *current = next;
  for (i = 0; i < nlocal; i++)
    current[i] = i;
  for (i = 0; i < nlocal; i++) {
    if (current[i] == permute[i])
      continue;
    avec->copy(i, nlocal, 0);
    empty = i;
    while (permute[empty] != i) {
      avec->copy(permute[empty], empty, 0);
      empty = current[empty] = permute[empty];
    }
    avec->copy(nlocal, empty, 0);
    current[empty] = permute[empty];
  }
}
void Atom::setup_sort_bins() {
  double binsize = 0.0;
  if (userbinsize > 0.0)
    binsize = userbinsize;
  else if (neighbor->cutneighmax > 0.0)
    binsize = 0.5 * neighbor->cutneighmax;
  if ((binsize == 0.0) && (sortfreq > 0)) {
    sortfreq = 0;
    if (comm->me == 0)
      error->warning(FLERR, "No pairwise cutoff or binsize set. Atom sorting "
                            "therefore disabled.");
    return;
  }
  double bininv = 1.0 / binsize;
  bboxlo[0] = domain->sublo[0];
  bboxlo[1] = domain->sublo[1];
  bboxlo[2] = domain->sublo[2];
  bboxhi[0] = domain->subhi[0];
  bboxhi[1] = domain->subhi[1];
  bboxhi[2] = domain->subhi[2];
  nbinx = static_cast<int>((bboxhi[0] - bboxlo[0]) * bininv);
  nbiny = static_cast<int>((bboxhi[1] - bboxlo[1]) * bininv);
  nbinz = static_cast<int>((bboxhi[2] - bboxlo[2]) * bininv);
  if (domain->dimension == 2)
    nbinz = 1;
  if (nbinx == 0)
    nbinx = 1;
  if (nbiny == 0)
    nbiny = 1;
  if (nbinz == 0)
    nbinz = 1;
  bininvx = nbinx / (bboxhi[0] - bboxlo[0]);
  bininvy = nbiny / (bboxhi[1] - bboxlo[1]);
  bininvz = nbinz / (bboxhi[2] - bboxlo[2]);
  if (1.0 * nbinx * nbiny * nbinz > INT_MAX)
    error->one(FLERR, "Too many atom sorting bins");
  nbins = nbinx * nbiny * nbinz;
  if (nbins > maxbin) {
    memory->destroy(binhead);
    maxbin = nbins;
    memory->create(binhead, maxbin, "atom:binhead");
  }
}
void Atom::update_callback(int ifix) {
  for (int i = 0; i < nextra_grow; i++)
    if (extra_grow[i] > ifix)
      extra_grow[i]--;
  for (int i = 0; i < nextra_restart; i++)
    if (extra_restart[i] > ifix)
      extra_restart[i]--;
  for (int i = 0; i < nextra_border; i++)
    if (extra_border[i] > ifix)
      extra_border[i]--;
}
