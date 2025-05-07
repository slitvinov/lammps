#include <set>
#include <map>
#include <vector>
#include <cstring>
#include <map>
#include <cstdio>
#include <string>
#include <mpi.h>
#include "lammps.h"
#include "pointers.h"
#include "force.h"
#include "lmptype.h"
#include "atom.h"
#include "comm.h"
#include "pair.h"
#include "pair_dpd.h"
#include "utils.h"
using namespace LAMMPS_NS;
template <typename S, typename T> static S *style_creator(LAMMPS *lmp) {
  return new T(lmp);
}
Force::Force(LAMMPS *lmp) : Pointers(lmp) {
  newton = newton_pair = 1;
  dielectric = 1.0;
  qqr2e_lammps_real = 332.06371;
  qqr2e_charmm_real = 332.0716;
  pair = nullptr;
  pair_style = utils::strdup("none");
  create_factories();
}
void _noopt Force::create_factories() {
  pair_map = new PairCreatorMap();
  (*pair_map)["dpd"] = &style_creator<Pair, PairDPD>;
}
Force::~Force() {
  delete[] pair_style;
  if (pair)
    delete pair;
  pair = nullptr;
  delete pair_map;
}
void Force::init() {
  qqrd2e = qqr2e / dielectric;
  if (pair)
    pair->init();
}
void Force::setup() {
  if (pair)
    pair->setup();
}
void Force::create_pair(const std::string &style, int trysuffix) {
  delete[] pair_style;
  if (pair)
    delete pair;
  pair_style = nullptr;
  pair = nullptr;
  int sflag;
  pair = new_pair(style, trysuffix, sflag);
  pair_style = store_style(style, sflag);
}
Pair *Force::new_pair(const std::string &style, int trysuffix, int &sflag) {
  sflag = 0;
  if (style == "none")
    return nullptr;
  if (pair_map->find(style) != pair_map->end()) {
    PairCreator &pair_creator = (*pair_map)[style];
    return pair_creator(lmp);
  }
  return nullptr;
}
Pair *Force::pair_match(const std::string &word, int exact, int nsub) {
  int iwhich, count;
  if (exact && (word == pair_style))
    return pair;
  else if (!exact && utils::strmatch(pair_style, word))
    return pair;
  return nullptr;
}
char *Force::pair_match_ptr(Pair *ptr) {
  if (ptr == pair)
    return pair_style;
  return nullptr;
}
char *Force::store_style(const std::string &style, int sflag) {
  std::string estyle = style;
  return utils::strdup(estyle);
}
