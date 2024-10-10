/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_FORCE_H
#define LMP_FORCE_H

#include "pointers.h"

#include <map>

namespace LAMMPS_NS {
class Angle;
class Pair;

enum { ENERGY_NONE = 0x00, ENERGY_GLOBAL = 0x01, ENERGY_ATOM = 0x02 };

// clang-format off
enum {
  VIRIAL_NONE     = 0x00,
  VIRIAL_PAIR     = 0x01,
  VIRIAL_FDOTR    = 0x02,
  VIRIAL_ATOM     = 0x04,
  VIRIAL_CENTROID = 0x08
};
// clang-format on

enum { CENTROID_SAME = 0, CENTROID_AVAIL = 1, CENTROID_NOTAVAIL = 2 };

class Force : protected Pointers {
 public:
  double boltz;          // Boltzmann constant (eng/degree-K)
  double hplanck;        // Planck's constant (energy-time)
  double mvv2e;          // conversion of mv^2 to energy
  double ftm2v;          // conversion of ft/m to velocity
  double mv2d;           // conversion of mass/volume to density
  double nktv2p;         // conversion of NkT/V to pressure
  double qqr2e;          // conversion of q^2/r to energy
  double qe2f;           // conversion of qE to force
  double vxmu2f;         // conversion of vx dynamic-visc to force
  double xxt2kmu;        // conversion of xx/t to kinematic-visc
  double dielectric;     // dielectric constant
  double qqrd2e;         // q^2/r to energy w/ dielectric constant
  double e_mass;         // electron mass
  double hhmrr2e;        // conversion of (hbar)^2/(mr^2) to energy
  double mvh2r;          // conversion of mv/hbar to distance
                         // hbar = h/(2*pi)
  double angstrom;       // 1 angstrom in native units
  double femtosecond;    // 1 femtosecond in native units
  double qelectron;      // 1 electron charge abs() in native units

  double qqr2e_lammps_real;    // different versions of this constant
  double qqr2e_charmm_real;    // used by new CHARMM pair styles

  int newton, newton_pair;    // Newton's 3rd law settings

  Pair *pair;
  char *pair_style;
  char *pair_restart;

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
  double memory_usage();

 private:
  void create_factories();
};

}    // namespace LAMMPS_NS

#endif
