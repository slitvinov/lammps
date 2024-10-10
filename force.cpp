/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "force.h"
#include "style_pair.h"        // IWYU pragma: keep
#include "atom.h"
#include "comm.h"
#include "error.h"

#include <cstring>

using namespace LAMMPS_NS;

// template for factory functions:
// there will be one instance for each style keyword in the respective style_xxx.h files

template <typename S, typename T> static S *style_creator(LAMMPS *lmp)
{
  return new T(lmp);
}

/* ---------------------------------------------------------------------- */

Force::Force(LAMMPS *lmp) : Pointers(lmp)
{
  newton = newton_pair = 1;

  dielectric = 1.0;
  qqr2e_lammps_real = 332.06371;    // these constants are toggled
  qqr2e_charmm_real = 332.0716;     // by new CHARMM pair styles

  pair = nullptr;
  pair_style = utils::strdup("none");
  pair_restart = nullptr;
  create_factories();
}

void _noopt Force::create_factories()
{
  // fill pair map with pair styles listed in style_pair.h

  pair_map = new PairCreatorMap();

#define PAIR_CLASS
#define PairStyle(key, Class) (*pair_map)[#key] = &style_creator<Pair, Class>;
#include "style_pair.h"    // IWYU pragma: keep
#undef PairStyle
#undef PAIR_CLASS
}

/* ---------------------------------------------------------------------- */

Force::~Force()
{
  delete[] pair_style;
  delete[] pair_restart;
  if (pair) delete pair;
  pair = nullptr;
  delete pair_map;
}

/* ---------------------------------------------------------------------- */

void Force::init()
{
  qqrd2e = qqr2e / dielectric;

  // check if pair style must be specified after restart
  if (pair_restart) {
    if (!pair)
      error->all(FLERR, "Must re-specify non-restarted pair style ({}) after read_restart",
                 pair_restart);
  }
  if (pair) pair->init();        // so g_ewald is defined
}

/* ---------------------------------------------------------------------- */

void Force::setup()
{
  if (pair) pair->setup();
}

/* ----------------------------------------------------------------------
   create a pair style, called from input script or restart file
------------------------------------------------------------------------- */

void Force::create_pair(const std::string &style, int trysuffix)
{
  delete[] pair_style;
  if (pair) delete pair;
  delete[] pair_restart;
  pair_style = nullptr;
  pair = nullptr;
  pair_restart = nullptr;

  int sflag;
  pair = new_pair(style, trysuffix, sflag);
  pair_style = store_style(style, sflag);
}

/* ----------------------------------------------------------------------
   generate a pair class
   if trysuffix = 1, try first with suffix1/2 appended
   return sflag = 0 for no suffix added, 1 or 2 for suffix1/2 added
------------------------------------------------------------------------- */

Pair *Force::new_pair(const std::string &style, int trysuffix, int &sflag)
{
  sflag = 0;
  if (style == "none") return nullptr;
  if (pair_map->find(style) != pair_map->end()) {
    PairCreator &pair_creator = (*pair_map)[style];
    return pair_creator(lmp);
  }
  return nullptr;
}

/* ----------------------------------------------------------------------
   return ptr to Pair class if matches word or matches hybrid sub-style
   if exact, then style name must be exact match to word
   if not exact, style name must contain word
   if nsub > 0, match Nth hybrid sub-style
   return nullptr if no match or if nsub=0 and multiple sub-styles match
------------------------------------------------------------------------- */

Pair *Force::pair_match(const std::string &word, int exact, int nsub)
{
  int iwhich, count;

  if (exact && (word == pair_style))
    return pair;
  else if (!exact && utils::strmatch(pair_style, word))
    return pair;
  return nullptr;
}

/* ----------------------------------------------------------------------
   return style name of Pair class that matches Pair ptr
   called by Neighbor::print_neigh_info()
   return nullptr if no match
------------------------------------------------------------------------- */

char *Force::pair_match_ptr(Pair *ptr)
{
  if (ptr == pair) return pair_style;
  return nullptr;
}


/* ----------------------------------------------------------------------
   store style name in str allocated here
   if sflag = 0, no suffix
   if sflag = 1/2, append suffix or suffix2 to style
------------------------------------------------------------------------- */

char *Force::store_style(const std::string &style, int sflag)
{
  std::string estyle = style;
  return utils::strdup(estyle);
}


/* ----------------------------------------------------------------------
   memory usage of force classes
------------------------------------------------------------------------- */

double Force::memory_usage()
{
  double bytes = 0;
  if (pair) bytes += pair->memory_usage();
  return bytes;
}
