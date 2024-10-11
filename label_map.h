#ifndef LMP_LABEL_MAP_H
#define LMP_LABEL_MAP_H 
#include "pointers.h"
#include <unordered_map>
namespace LAMMPS_NS {
class LabelMap : protected Pointers {
  friend class AtomVec;
  friend class ReadData;
 public:
  LabelMap(LAMMPS *lmp, int, int, int, int, int);
  ~LabelMap();
  void modify_lmap(int, char **);
  void merge_lmap(LabelMap *, int);
  void create_lmap2lmap(LabelMap *, int);
  int find(const std::string &, int) const;
  bool is_complete(int) const;
  void write_data(FILE *);
  void read_restart(FILE *fp);
  void write_restart(FILE *);
 protected:
  int natomtypes, nbondtypes, nangletypes, ndihedraltypes, nimpropertypes;
  std::vector<std::string> typelabel, btypelabel, atypelabel;
  std::vector<std::string> dtypelabel, itypelabel;
  std::unordered_map<std::string, int> typelabel_map;
  std::unordered_map<std::string, int> btypelabel_map;
  std::unordered_map<std::string, int> atypelabel_map;
  std::unordered_map<std::string, int> dtypelabel_map;
  std::unordered_map<std::string, int> itypelabel_map;
  struct Lmap2Lmap {
    int *atom;
    int *bond;
    int *angle;
    int *dihedral;
    int *improper;
  };
  Lmap2Lmap lmap2lmap;
  void reset_type_labels();
  int find_or_create(const std::string &, std::vector<std::string> &,
                     std::unordered_map<std::string, int> &);
  int search(const std::string &,
             const std::unordered_map<std::string, int> &) const;
  char *read_string(FILE *);
  void write_string(const std::string &, FILE *);
  int read_int(FILE *);
  void write_map(const std::string &);
};
}
#endif
