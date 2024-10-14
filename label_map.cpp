#include "label_map.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include <cstring>
using namespace LAMMPS_NS;
LabelMap::LabelMap(LAMMPS *_lmp, int _natomtypes, int _nbondtypes,
                   int _nangletypes, int _ndihedraltypes, int _nimpropertypes)
    : Pointers(_lmp), natomtypes(_natomtypes), nbondtypes(_nbondtypes),
      nangletypes(_nangletypes), ndihedraltypes(_ndihedraltypes),
      nimpropertypes(_nimpropertypes) {
  lmap2lmap.atom = lmap2lmap.bond = lmap2lmap.angle = lmap2lmap.dihedral =
      lmap2lmap.improper = nullptr;
  reset_type_labels();
}
LabelMap::~LabelMap() {
  delete[] lmap2lmap.atom;
  delete[] lmap2lmap.bond;
  delete[] lmap2lmap.angle;
  delete[] lmap2lmap.dihedral;
  delete[] lmap2lmap.improper;
}
void LabelMap::reset_type_labels() {
  typelabel_map.clear();
  typelabel.resize(natomtypes);
  delete[] lmap2lmap.atom;
  lmap2lmap.atom = new int[natomtypes];
  for (auto &i : typelabel)
    i.clear();
  memset(lmap2lmap.atom, 0, natomtypes * sizeof(int));
  btypelabel_map.clear();
  btypelabel.resize(nbondtypes);
  delete[] lmap2lmap.bond;
  for (auto &i : btypelabel)
    i.clear();
  lmap2lmap.bond = new int[nbondtypes];
  memset(lmap2lmap.bond, 0, nbondtypes * sizeof(int));
  atypelabel_map.clear();
  atypelabel.resize(nangletypes);
  delete[] lmap2lmap.angle;
  for (auto &i : atypelabel)
    i.clear();
  lmap2lmap.angle = new int[nangletypes];
  memset(lmap2lmap.angle, 0, nangletypes * sizeof(int));
  dtypelabel_map.clear();
  dtypelabel.resize(ndihedraltypes);
  delete[] lmap2lmap.dihedral;
  for (auto &i : dtypelabel)
    i.clear();
  lmap2lmap.dihedral = new int[ndihedraltypes];
  memset(lmap2lmap.dihedral, 0, ndihedraltypes * sizeof(int));
  itypelabel_map.clear();
  itypelabel.resize(nimpropertypes);
  delete[] lmap2lmap.improper;
  for (auto &i : itypelabel)
    i.clear();
  lmap2lmap.improper = new int[nimpropertypes];
  memset(lmap2lmap.improper, 0, nimpropertypes * sizeof(int));
}
void LabelMap::modify_lmap(int narg, char **arg) {
  if ((narg < 1) || ((narg > 2) && ((narg % 2) == 0)))
    error->all(FLERR, "Incorrect number of arguments for labelmap command");
  int ntypes;
  std::vector<std::string> *labels;
  std::unordered_map<std::string, int> *labels_map;
  const std::string tlabel(arg[0]);
  if (tlabel == "atom") {
    ntypes = natomtypes;
    labels = &typelabel;
    labels_map = &typelabel_map;
  } else if (tlabel == "bond") {
    ntypes = nbondtypes;
    labels = &btypelabel;
    labels_map = &btypelabel_map;
  } else if (tlabel == "angle") {
    ntypes = nangletypes;
    labels = &atypelabel;
    labels_map = &atypelabel_map;
  } else if (tlabel == "dihedral") {
    ntypes = ndihedraltypes;
    labels = &dtypelabel;
    labels_map = &dtypelabel_map;
  } else if (tlabel == "improper") {
    ntypes = nimpropertypes;
    labels = &itypelabel;
    labels_map = &itypelabel_map;
  } else if (tlabel == "clear") {
    if (narg != 1)
      error->all(FLERR,
                 "Incorrect number of arguments for labelmap clear command");
    reset_type_labels();
    return;
  } else if (tlabel == "write") {
    if (narg != 2)
      error->all(FLERR,
                 "Incorrect number of arguments for labelmap write command");
    write_map(arg[1]);
    return;
  } else
    error->all(FLERR, "Unknown labelmap keyword {}", tlabel);
  int iarg = 1;
  if (narg == 1)
    utils::missing_cmd_args(FLERR, "labelmap " + tlabel, error);
  while (iarg < narg) {
    if (iarg + 2 > narg)
      utils::missing_cmd_args(FLERR, "labelmap " + tlabel, error);
    if (ntypes < 1)
      error->all(FLERR, "No {} types allowed with current box settings",
                 tlabel);
    int itype = utils::inumeric(FLERR, arg[iarg++], false, lmp);
    if ((itype < 1) || (itype > ntypes))
      error->all(FLERR, "Labelmap {} type {} must be within 1-{}", tlabel,
                 itype, ntypes);
    std::string slabel = utils::utf8_subst(utils::trim(arg[iarg++]));
    if (utils::is_type(slabel) != 1)
      error->all(FLERR, "Type label string {} for {} type {} is invalid",
                 slabel, tlabel, itype);
    int found = search(slabel, (*labels_map));
    if ((found != -1) && (found != itype))
      error->all(FLERR, "The {} type label {} is already in use for type {}",
                 tlabel, slabel, (*labels_map)[slabel]);
    std::string &str = (*labels)[itype - 1];
    if (!str.empty())
      (*labels_map).erase(str);
    str = slabel;
    (*labels_map)[slabel] = itype;
  }
}
void LabelMap::merge_lmap(LabelMap *lmap2, int mode) {
  switch (mode) {
  case Atom::ATOM:
    for (auto &it : lmap2->typelabel)
      find_or_create(it, typelabel, typelabel_map);
    break;
  }
}
void LabelMap::create_lmap2lmap(LabelMap *lmap2, int mode) {
  switch (mode) {
  case Atom::ATOM:
    for (int i = 0; i < natomtypes; ++i)
      lmap2lmap.atom[i] = search(typelabel[i], lmap2->typelabel_map);
    break;
  }
}
int LabelMap::find_or_create(const std::string &mylabel,
                             std::vector<std::string> &labels,
                             std::unordered_map<std::string, int> &labels_map) {
  auto search = labels_map.find(mylabel);
  if (search != labels_map.end())
    return search->second;
  auto labels_map_size = labels_map.size();
  if (labels_map_size < labels.size()) {
    labels[labels_map_size] = mylabel;
    int index = static_cast<int>(labels_map_size + 1);
    labels_map[mylabel] = index;
    return index;
  }
  error->all(FLERR, "Topology type exceeds system topology type");
  return -1;
}
int LabelMap::find(const std::string &mylabel, int mode) const {
  switch (mode) {
  case Atom::ATOM:
    return search(mylabel, typelabel_map);
    break;
  default:
    return -1;
  }
}
int LabelMap::search(
    const std::string &mylabel,
    const std::unordered_map<std::string, int> &labels_map) const {
  auto search = labels_map.find(mylabel);
  if (search == labels_map.end())
    return -1;
  return search->second;
}
bool LabelMap::is_complete(int mode) const {
  switch (mode) {
  case Atom::ATOM:
    return static_cast<int>(typelabel_map.size()) == natomtypes;
    break;
  }
  return false;
}
char *LabelMap::read_string(FILE *fp) {
  int n = read_int(fp);
  if (n < 0)
    error->all(FLERR, "Illegal size string or corrupt restart");
  char *value = new char[n];
  if (comm->me == 0)
    utils::sfread(FLERR, value, sizeof(char), n, fp, nullptr, error);
  MPI_Bcast(value, n, MPI_CHAR, 0, world);
  return value;
}
void LabelMap::write_string(const std::string &str, FILE *fp) {
  const char *cstr = str.c_str();
  int n = strlen(cstr) + 1;
  fwrite(&n, sizeof(int), 1, fp);
  fwrite(cstr, sizeof(char), n, fp);
}
int LabelMap::read_int(FILE *fp) {
  int value;
  if ((comm->me == 0) && (fread(&value, sizeof(int), 1, fp) < 1))
    value = -1;
  MPI_Bcast(&value, 1, MPI_INT, 0, world);
  return value;
}
void LabelMap::write_map(const std::string &filename) {
  if (comm->me == 0) {
    FILE *fp = fopen(filename.c_str(), "w");
    if (!fp)
      error->one(FLERR, "Cannot open label map file {}: {}", filename,
                 utils::getsyserror());
    if (typelabel_map.size() > 0) {
      fputs("labelmap atom", fp);
      for (int i = 0; i < natomtypes; ++i)
        if (!typelabel[i].empty())
          fmt::print(fp, " {} \"\"\" {} \"\"\"", i + 1, typelabel[i]);
      fputc('\n', fp);
    }
    if (btypelabel_map.size() > 0) {
      fputs("labelmap bond", fp);
      for (int i = 0; i < nbondtypes; ++i)
        if (!btypelabel[i].empty())
          fmt::print(fp, " {} \"\"\" {} \"\"\"", i + 1, btypelabel[i]);
      fputc('\n', fp);
    }
    if (atypelabel_map.size() > 0) {
      fputs("labelmap angle", fp);
      for (int i = 0; i < nangletypes; ++i)
        if (!atypelabel[i].empty())
          fmt::print(fp, " {} \"\"\" {} \"\"\"", i + 1, atypelabel[i]);
      fputc('\n', fp);
    }
    if (dtypelabel_map.size() > 0) {
      fputs("labelmap dihedral", fp);
      for (int i = 0; i < ndihedraltypes; ++i)
        if (!dtypelabel[i].empty())
          fmt::print(fp, " {} \"\"\" {} \"\"\"", i + 1, dtypelabel[i]);
      fputc('\n', fp);
    }
    if (itypelabel_map.size() > 0) {
      fputs("labelmap improper", fp);
      for (int i = 0; i < nimpropertypes; ++i)
        if (!itypelabel[i].empty())
          fmt::print(fp, " {} \"\"\" {} \"\"\"", i + 1, itypelabel[i]);
      fputc('\n', fp);
    }
    fclose(fp);
  }
}
