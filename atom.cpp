#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "label_map.h"
#include "library.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "style_atom.h"
#include "tokenizer.h"
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
void Atom::add_peratom_vary(const std::string &name, void *address,
                            int datatype, int *cols, void *length,
                            int collength) {
  PerAtom item = {name, address, length, cols, datatype, -1, collength, 0};
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
AtomVec *Atom::style_match(const char *style) {
  if (strcmp(atom_style, style) == 0)
    return avec;
  return nullptr;
}
void Atom::modify_params(int narg, char **arg) {
  if (narg == 0)
    utils::missing_cmd_args(FLERR, "atom_modify", error);
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "id") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, "atom_modify id", error);
      if (domain->box_exist)
        error->all(FLERR,
                   "Atom_modify id command after simulation box is defined");
      tag_enable = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "map") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, "atom_modify map", error);
      if (domain->box_exist)
        error->all(FLERR,
                   "Atom_modify map command after simulation box is defined");
      if (strcmp(arg[iarg + 1], "array") == 0)
        map_user = 1;
      else if (strcmp(arg[iarg + 1], "hash") == 0)
        map_user = 2;
      else if (strcmp(arg[iarg + 1], "yes") == 0)
        map_user = 3;
      else
        error->all(FLERR, "Illegal atom_modify map command argument {}",
                   arg[iarg + 1]);
      map_style = map_user;
      iarg += 2;
    } else if (strcmp(arg[iarg], "first") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, "atom_modify first", error);
      if (strcmp(arg[iarg + 1], "all") == 0) {
        delete[] firstgroupname;
        firstgroupname = nullptr;
      } else {
        firstgroupname = utils::strdup(arg[iarg + 1]);
        sortfreq = 0;
      }
      iarg += 2;
    } else if (strcmp(arg[iarg], "sort") == 0) {
      if (iarg + 3 > narg)
        utils::missing_cmd_args(FLERR, "atom_modify sort", error);
      sortfreq = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      userbinsize = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      if (sortfreq < 0)
        error->all(FLERR, "Illegal atom_modify sort frequency {}", sortfreq);
      if (userbinsize < 0.0)
        error->all(FLERR, "Illegal atom_modify sort bin size {}", userbinsize);
      if ((sortfreq >= 0) && firstgroupname)
        error->all(
            FLERR,
            "Atom_modify sort and first options cannot be used together");
      iarg += 3;
    } else
      error->all(FLERR, "Illegal atom_modify command argument: {}", arg[iarg]);
  }
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
void Atom::set_mass(const char *file, int line, const char *str,
                    int type_offset, int labelflag, int *ilabel) {
  if (mass == nullptr)
    error->all(file, line, "Cannot set mass for atom style {}", atom_style);
  int itype;
  double mass_one;
  auto location = "Masses section of data file";
  auto values = Tokenizer(str).as_vector();
  int nwords = values.size();
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (utils::strmatch(values[i], "^#")) {
      nwords = i;
      break;
    }
  }
  if (nwords != 2)
    error->all(file, line, "Invalid format in {}: {}", location, str);
  auto typestr = utils::utf8_subst(values[0]);
  switch (utils::is_type(typestr)) {
  case 0: {
    itype = utils::inumeric(file, line, typestr, false, lmp);
    if ((itype < 1) || (itype > ntypes))
      error->all(file, line, "Invalid atom type {} in {}: {}", itype, location,
                 utils::trim(str));
    if (labelflag)
      itype = ilabel[itype - 1];
    break;
  }
  case 1: {
    error->all(file, line, "Invalid atom type in {}: {}", location,
               utils::trim(str));
    break;
  }
  default:
    itype = -1000000000;
    error->one(file, line, "Invalid {}: {}", location, utils::trim(str));
    break;
  }
  itype += type_offset;
  mass_one = utils::numeric(file, line, values[1], false, lmp);
  if (itype < 1 || itype > ntypes)
    error->all(file, line, "Invalid atom type {} in {}: {}", itype, location,
               utils::trim(str));
  if (mass_one <= 0.0)
    error->all(file, line, "Invalid mass value {} in {}: {}", mass_one,
               location, utils::trim(str));
  mass[itype] = mass_one;
  mass_setflag[itype] = 1;
}
void Atom::set_mass(const char *file, int line, int itype, double value) {
  if (mass == nullptr)
    error->all(file, line, "Cannot set per-type mass for atom style {}",
               atom_style);
  if (itype < 1 || itype > ntypes)
    error->all(file, line, "Invalid type {} for atom mass {}", itype, value);
  if (value <= 0.0) {
    if (comm->me == 0)
      error->warning(file, line,
                     "Ignoring invalid mass value {} for atom type {}", value,
                     itype);
    return;
  }
  mass[itype] = value;
  mass_setflag[itype] = 1;
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
void Atom::set_mass(double *values) {
  for (int itype = 1; itype <= ntypes; itype++) {
    mass[itype] = values[itype];
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
int Atom::radius_consistency(int itype, double &rad) {
  double value = -1.0;
  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    if (type[i] != itype)
      continue;
    if (value < 0.0)
      value = radius[i];
    else if (value != radius[i])
      flag = 1;
  }
  int flagall;
  MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_SUM, world);
  if (flagall)
    return 0;
  MPI_Allreduce(&value, &rad, 1, MPI_DOUBLE, MPI_MAX, world);
  return 1;
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
  if (domain->triclinic)
    domain->bbox(domain->sublo_lamda, domain->subhi_lamda, bboxlo, bboxhi);
  else {
    bboxlo[0] = domain->sublo[0];
    bboxlo[1] = domain->sublo[1];
    bboxlo[2] = domain->sublo[2];
    bboxhi[0] = domain->subhi[0];
    bboxhi[1] = domain->subhi[1];
    bboxhi[2] = domain->subhi[2];
  }
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
void Atom::add_callback(int flag) {
  int ifix;
  for (ifix = 0; ifix < modify->nfix; ifix++)
    if (modify->fix[ifix] == nullptr)
      break;
  if (flag == GROW) {
    if (nextra_grow == nextra_grow_max) {
      nextra_grow_max += DELTA;
      memory->grow(extra_grow, nextra_grow_max, "atom:extra_grow");
    }
    extra_grow[nextra_grow] = ifix;
    nextra_grow++;
    std::sort(extra_grow, extra_grow + nextra_grow);
  } else if (flag == RESTART) {
    if (nextra_restart == nextra_restart_max) {
      nextra_restart_max += DELTA;
      memory->grow(extra_restart, nextra_restart_max, "atom:extra_restart");
    }
    extra_restart[nextra_restart] = ifix;
    nextra_restart++;
    std::sort(extra_restart, extra_restart + nextra_restart);
  } else if (flag == BORDER) {
    if (nextra_border == nextra_border_max) {
      nextra_border_max += DELTA;
      memory->grow(extra_border, nextra_border_max, "atom:extra_border");
    }
    extra_border[nextra_border] = ifix;
    nextra_border++;
    std::sort(extra_border, extra_border + nextra_border);
  }
}
void Atom::delete_callback(const char *id, int flag) {
  if (id == nullptr)
    return;
  int ifix = modify->find_fix(id);
  if (flag == GROW) {
    int match;
    for (match = 0; match < nextra_grow; match++)
      if (extra_grow[match] == ifix)
        break;
    if ((nextra_grow == 0) || (match == nextra_grow))
      error->all(FLERR, "Trying to delete non-existent Atom::grow() callback");
    for (int i = match; i < nextra_grow - 1; i++)
      extra_grow[i] = extra_grow[i + 1];
    nextra_grow--;
  } else if (flag == RESTART) {
    int match;
    for (match = 0; match < nextra_restart; match++)
      if (extra_restart[match] == ifix)
        break;
    if ((nextra_restart == 0) || (match == nextra_restart))
      error->all(FLERR,
                 "Trying to delete non-existent Atom::restart() callback");
    for (int i = match; i < nextra_restart - 1; i++)
      extra_restart[i] = extra_restart[i + 1];
    nextra_restart--;
  } else if (flag == BORDER) {
    int match;
    for (match = 0; match < nextra_border; match++)
      if (extra_border[match] == ifix)
        break;
    if ((nextra_border == 0) || (match == nextra_border))
      error->all(FLERR,
                 "Trying to delete non-existent Atom::border() callback");
    for (int i = match; i < nextra_border - 1; i++)
      extra_border[i] = extra_border[i + 1];
    nextra_border--;
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
int Atom::find_custom(const char *name, int &flag, int &cols) {
  if (name == nullptr)
    return -1;
  return -1;
}
int Atom::add_custom(const char *name, int flag, int cols) {
  int index = -1;
  if ((flag == 0) && (cols == 0)) {
    index = nivector;
    nivector++;
    ivname = (char **)memory->srealloc(ivname, nivector * sizeof(char *),
                                       "atom:ivname");
    ivname[index] = utils::strdup(name);
    ivector = (int **)memory->srealloc(ivector, nivector * sizeof(int *),
                                       "atom:ivector");
    memory->create(ivector[index], nmax, "atom:ivector");
  } else if ((flag == 1) && (cols == 0)) {
    index = ndvector;
    ndvector++;
    dvname = (char **)memory->srealloc(dvname, ndvector * sizeof(char *),
                                       "atom:dvname");
    dvname[index] = utils::strdup(name);
    dvector = (double **)memory->srealloc(dvector, ndvector * sizeof(double *),
                                          "atom:dvector");
    memory->create(dvector[index], nmax, "atom:dvector");
  } else if ((flag == 0) && (cols > 0)) {
    index = niarray;
    niarray++;
    ianame = (char **)memory->srealloc(ianame, niarray * sizeof(char *),
                                       "atom:ianame");
    ianame[index] = utils::strdup(name);
    iarray = (int ***)memory->srealloc(iarray, niarray * sizeof(int **),
                                       "atom:iarray");
    memory->create(iarray[index], nmax, cols, "atom:iarray");
    icols = (int *)memory->srealloc(icols, niarray * sizeof(int), "atom:icols");
    icols[index] = cols;
  } else if ((flag == 1) && (cols > 0)) {
    index = ndarray;
    ndarray++;
    daname = (char **)memory->srealloc(daname, ndarray * sizeof(char *),
                                       "atom:daname");
    daname[index] = utils::strdup(name);
    darray = (double ***)memory->srealloc(darray, ndarray * sizeof(double **),
                                          "atom:darray");
    memory->create(darray[index], nmax, cols, "atom:darray");
    dcols = (int *)memory->srealloc(dcols, ndarray * sizeof(int), "atom:dcols");
    dcols[index] = cols;
  }
  if (index < 0)
    error->all(FLERR, "Invalid call to Atom::add_custom()");
  return index;
}
void Atom::remove_custom(int index, int flag, int cols) {
  if (flag == 0 && cols == 0) {
    memory->destroy(ivector[index]);
    ivector[index] = nullptr;
    delete[] ivname[index];
    ivname[index] = nullptr;
  } else if (flag == 1 && cols == 0) {
    memory->destroy(dvector[index]);
    dvector[index] = nullptr;
    delete[] dvname[index];
    dvname[index] = nullptr;
  } else if (flag == 0 && cols) {
    memory->destroy(iarray[index]);
    iarray[index] = nullptr;
    delete[] ianame[index];
    ianame[index] = nullptr;
  } else if (flag == 1 && cols) {
    memory->destroy(darray[index]);
    darray[index] = nullptr;
    delete[] daname[index];
    daname[index] = nullptr;
  }
}
void *Atom::extract(const char *name) {
  if (strcmp(name, "mass") == 0)
    return (void *)mass;
  if (strcmp(name, "id") == 0)
    return (void *)tag;
  if (strcmp(name, "type") == 0)
    return (void *)type;
  if (strcmp(name, "mask") == 0)
    return (void *)mask;
  if (strcmp(name, "image") == 0)
    return (void *)image;
  if (strcmp(name, "x") == 0)
    return (void *)x;
  if (strcmp(name, "v") == 0)
    return (void *)v;
  if (strcmp(name, "f") == 0)
    return (void *)f;
  if (strcmp(name, "q") == 0)
    return (void *)q;
  if (strcmp(name, "mu") == 0)
    return (void *)mu;
  if (strcmp(name, "omega") == 0)
    return (void *)omega;
  if (strcmp(name, "angmom") == 0)
    return (void *)angmom;
  if (strcmp(name, "torque") == 0)
    return (void *)torque;
  if (strcmp(name, "radius") == 0)
    return (void *)radius;
  if (strcmp(name, "rmass") == 0)
    return (void *)rmass;
  if (strcmp(name, "quat") == 0)
    return (void *)quat;
  if (strcmp(name, "temperature") == 0)
    return (void *)temperature;
  if (strcmp(name, "heatflow") == 0)
    return (void *)heatflow;
  if (strcmp(name, "vfrac") == 0)
    return (void *)vfrac;
  if (strcmp(name, "s0") == 0)
    return (void *)s0;
  if (strcmp(name, "x0") == 0)
    return (void *)x0;
  if (strcmp(name, "sp") == 0)
    return (void *)sp;
  if (strcmp(name, "espin") == 0)
    return (void *)spin;
  if (strcmp(name, "spin") == 0)
    return (void *)spin;
  if (strcmp(name, "eradius") == 0)
    return (void *)eradius;
  if (strcmp(name, "ervel") == 0)
    return (void *)ervel;
  if (strcmp(name, "erforce") == 0)
    return (void *)erforce;
  if (strcmp(name, "ervelforce") == 0)
    return (void *)ervelforce;
  if (strcmp(name, "cs") == 0)
    return (void *)cs;
  if (strcmp(name, "csforce") == 0)
    return (void *)csforce;
  if (strcmp(name, "vforce") == 0)
    return (void *)vforce;
  if (strcmp(name, "etag") == 0)
    return (void *)etag;
  if (strcmp(name, "rho") == 0)
    return (void *)rho;
  if (strcmp(name, "drho") == 0)
    return (void *)drho;
  if (strcmp(name, "esph") == 0)
    return (void *)esph;
  if (strcmp(name, "desph") == 0)
    return (void *)desph;
  if (strcmp(name, "cv") == 0)
    return (void *)cv;
  if (strcmp(name, "vest") == 0)
    return (void *)vest;
  if (strcmp(name, "contact_radius") == 0)
    return (void *)contact_radius;
  if (strcmp(name, "smd_data_9") == 0)
    return (void *)smd_data_9;
  if (strcmp(name, "smd_stress") == 0)
    return (void *)smd_stress;
  if (strcmp(name, "eff_plastic_strain") == 0)
    return (void *)eff_plastic_strain;
  if (strcmp(name, "eff_plastic_strain_rate") == 0)
    return (void *)eff_plastic_strain_rate;
  if (strcmp(name, "damage") == 0)
    return (void *)damage;
  if (strcmp(name, "edpd_temp") == 0)
    return (void *)edpd_temp;
  if (utils::strmatch(name, "^[id]2?_")) {
    int which = 0, array = 0;
    if (name[0] == 'd')
      which = 1;
    if (name[1] == '2')
      array = 1;
    int index, flag, cols;
    if (!array)
      index = find_custom(&name[2], flag, cols);
    else
      index = find_custom(&name[3], flag, cols);
    if (index < 0)
      return nullptr;
    if (which != flag)
      return nullptr;
    if ((!array && cols) || (array && !cols))
      return nullptr;
    if (!which && !array)
      return (void *)ivector[index];
    if (which && !array)
      return (void *)dvector[index];
    if (!which && array)
      return (void *)iarray[index];
    if (which && array)
      return (void *)darray[index];
  }
  return nullptr;
}
int Atom::extract_datatype(const char *name) {
  if (strcmp(name, "mass") == 0)
    return LAMMPS_DOUBLE;
  if (strcmp(name, "id") == 0)
    return LAMMPS_TAGINT;
  if (strcmp(name, "type") == 0)
    return LAMMPS_INT;
  if (strcmp(name, "mask") == 0)
    return LAMMPS_INT;
  if (strcmp(name, "image") == 0)
    return LAMMPS_TAGINT;
  if (strcmp(name, "x") == 0)
    return LAMMPS_DOUBLE_2D;
  if (strcmp(name, "v") == 0)
    return LAMMPS_DOUBLE_2D;
  if (strcmp(name, "f") == 0)
    return LAMMPS_DOUBLE_2D;
  if (strcmp(name, "molecule") == 0)
    return LAMMPS_TAGINT;
  if (strcmp(name, "q") == 0)
    return LAMMPS_DOUBLE;
  if (strcmp(name, "mu") == 0)
    return LAMMPS_DOUBLE_2D;
  if (strcmp(name, "omega") == 0)
    return LAMMPS_DOUBLE_2D;
  if (strcmp(name, "angmom") == 0)
    return LAMMPS_DOUBLE_2D;
  if (strcmp(name, "torque") == 0)
    return LAMMPS_DOUBLE_2D;
  if (strcmp(name, "radius") == 0)
    return LAMMPS_DOUBLE;
  if (strcmp(name, "rmass") == 0)
    return LAMMPS_DOUBLE;
  if (strcmp(name, "ellipsoid") == 0)
    return LAMMPS_INT;
  if (strcmp(name, "line") == 0)
    return LAMMPS_INT;
  if (strcmp(name, "tri") == 0)
    return LAMMPS_INT;
  if (strcmp(name, "body") == 0)
    return LAMMPS_INT;
  if (strcmp(name, "quat") == 0)
    return LAMMPS_DOUBLE_2D;
  if (strcmp(name, "temperature") == 0)
    return LAMMPS_DOUBLE;
  if (strcmp(name, "heatflow") == 0)
    return LAMMPS_DOUBLE;
  if (strcmp(name, "vfrac") == 0)
    return LAMMPS_DOUBLE;
  if (strcmp(name, "s0") == 0)
    return LAMMPS_DOUBLE;
  if (strcmp(name, "x0") == 0)
    return LAMMPS_DOUBLE_2D;
  if (strcmp(name, "espin") == 0)
    return LAMMPS_INT;
  if (strcmp(name, "spin") == 0)
    return LAMMPS_INT;
  if (strcmp(name, "eradius") == 0)
    return LAMMPS_DOUBLE;
  if (strcmp(name, "ervel") == 0)
    return LAMMPS_DOUBLE;
  if (strcmp(name, "erforce") == 0)
    return LAMMPS_DOUBLE;
  if (strcmp(name, "ervelforce") == 0)
    return LAMMPS_DOUBLE;
  if (strcmp(name, "cs") == 0)
    return LAMMPS_DOUBLE_2D;
  if (strcmp(name, "csforce") == 0)
    return LAMMPS_DOUBLE_2D;
  if (strcmp(name, "vforce") == 0)
    return LAMMPS_DOUBLE_2D;
  if (strcmp(name, "etag") == 0)
    return LAMMPS_INT;
  if (strcmp(name, "rho") == 0)
    return LAMMPS_DOUBLE;
  if (strcmp(name, "drho") == 0)
    return LAMMPS_DOUBLE;
  if (strcmp(name, "esph") == 0)
    return LAMMPS_DOUBLE;
  if (strcmp(name, "desph") == 0)
    return LAMMPS_DOUBLE;
  if (strcmp(name, "cv") == 0)
    return LAMMPS_DOUBLE;
  if (strcmp(name, "vest") == 0)
    return LAMMPS_DOUBLE_2D;
  if (strcmp(name, "contact_radius") == 0)
    return LAMMPS_DOUBLE;
  if (strcmp(name, "smd_data_9") == 0)
    return LAMMPS_DOUBLE_2D;
  if (strcmp(name, "smd_stress") == 0)
    return LAMMPS_DOUBLE_2D;
  if (strcmp(name, "eff_plastic_strain") == 0)
    return LAMMPS_DOUBLE;
  if (strcmp(name, "eff_plastic_strain_rate") == 0)
    return LAMMPS_DOUBLE;
  if (strcmp(name, "damage") == 0)
    return LAMMPS_DOUBLE;
  if (strcmp(name, "dpdTheta") == 0)
    return LAMMPS_DOUBLE;
  if (strcmp(name, "edpd_temp") == 0)
    return LAMMPS_DOUBLE;
  if (strcmp(name, "area") == 0)
    return LAMMPS_DOUBLE;
  if (strcmp(name, "ed") == 0)
    return LAMMPS_DOUBLE;
  if (strcmp(name, "em") == 0)
    return LAMMPS_DOUBLE;
  if (strcmp(name, "epsilon") == 0)
    return LAMMPS_DOUBLE;
  if (strcmp(name, "curvature") == 0)
    return LAMMPS_DOUBLE;
  if (strcmp(name, "q_unscaled") == 0)
    return LAMMPS_DOUBLE;
  if (utils::strmatch(name, "^[id]2?_")) {
    int which = 0, array = 0;
    if (name[0] == 'd')
      which = 1;
    if (name[1] == '2')
      array = 1;
    int index, flag, cols;
    if (!array)
      index = find_custom(&name[2], flag, cols);
    else
      index = find_custom(&name[3], flag, cols);
    if (index < 0)
      return -1;
    if (which != flag)
      return -1;
    if ((!array && cols) || (array && !cols))
      return -1;
    if (which == 0)
      return LAMMPS_INT;
    else
      return LAMMPS_DOUBLE;
  }
  return -1;
}
