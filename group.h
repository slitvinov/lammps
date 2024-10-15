#ifndef LMP_GROUP_H
#define LMP_GROUP_H
#include "pointers.h"
#include <map>
namespace LAMMPS_NS {
class Region;
class Group : protected Pointers {
public:
  int ngroup;
  char **names;
  int *bitmask;
  int *inversemask;
  int *dynamic;
  Group(class LAMMPS *);
  ~Group() override;
  int find(const std::string &);

private:
  int me;
  std::map<tagint, int> *hash;
  int molbit;
};
} // namespace LAMMPS_NS
#endif
