#ifndef LMP_FORCE_H
#define LMP_FORCE_H
#include "pointers.h"
#include <map>
namespace LAMMPS_NS {
class Angle;
class Pair;
enum { ENERGY_NONE = 0x00, ENERGY_GLOBAL = 0x01, ENERGY_ATOM = 0x02 };
enum {
  VIRIAL_NONE = 0x00,
  VIRIAL_PAIR = 0x01,
  VIRIAL_FDOTR = 0x02,
  VIRIAL_ATOM = 0x04,
  VIRIAL_CENTROID = 0x08
};
enum { CENTROID_SAME = 0, CENTROID_AVAIL = 1, CENTROID_NOTAVAIL = 2 };
class Force : protected Pointers {
public:
  double boltz;
  double hplanck;
  double mvv2e;
  double ftm2v;
  double mv2d;
  double nktv2p;
  double qqr2e;
  double qe2f;
  double vxmu2f;
  double xxt2kmu;
  double dielectric;
  double qqrd2e;
  double e_mass;
  double hhmrr2e;
  double mvh2r;
  double angstrom;
  double femtosecond;
  double qelectron;
  double qqr2e_lammps_real;
  double qqr2e_charmm_real;
  int newton, newton_pair;
  Pair *pair;
  char *pair_style;
  typedef Pair *(*PairCreator)(LAMMPS *);
  typedef std::map<std::string, PairCreator> PairCreatorMap;
  PairCreatorMap *pair_map;
  Force(class LAMMPS *);
  ~Force() override;
  void init();
  void setup();
  void create_pair(const std::string &, int);
  Pair *new_pair(const std::string &, int, int &);
  Pair *pair_match(const std::string &, int, int nsub = 0);
  char *pair_match_ptr(Pair *);
  char *store_style(const std::string &, int);

private:
  void create_factories();
};
} // namespace LAMMPS_NS
#endif
