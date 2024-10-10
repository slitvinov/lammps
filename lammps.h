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

#ifndef LMP_LAMMPS_H
#define LMP_LAMMPS_H

#include <cstdio>
#include <mpi.h>

namespace LAMMPS_NS {

class LAMMPS {
 public:
  // ptrs to fundamental LAMMPS classes
  class Memory *memory;            // memory allocation functions
  class Error *error;              // error handling
  class Universe *universe;        // universe of processors
  class Input *input;              // input script processing
                                   // ptrs to top-level LAMMPS-specific classes
  class Atom *atom;                // atom-based quantities
  class Update *update;            // integrators/minimizers
  class Neighbor *neighbor;        // neighbor lists
  class Comm *comm;                // inter-processor communication
  class Domain *domain;            // simulation box
  class Force *force;              // inter-particle forces
  class Modify *modify;            // fixes and computes
  class Group *group;              // groups of atoms
                                   //
  const char *version;    // LAMMPS version string = date
  int num_ver;            // numeric version id derived from *version*
                          // that is constructed so that will be greater
                          // for newer versions in numeric or string
                          // value comparisons
  int restart_ver;        // -1 or numeric version id of LAMMPS version in restart
                          // file, in case LAMMPS was initialized from a restart
                          //
  MPI_Comm world;         // MPI communicator
  FILE *infile;           // infile
  FILE *screen;           // screen output
  FILE *logfile;          // logfile
                          //
  double initclock;       // wall clock at instantiation
  int skiprunflag;        // 1 inserts timer command to skip run and minimize loops
  int pair_only_flag;        // 1 if only force field pair styles are accelerated, 0 if all
  char *exename;             // pointer to argv[0]

  char ***packargs;    // arguments for cmdline package commands
  int num_package;     // number of cmdline package commands

  MPI_Comm external_comm;    // MPI comm encompassing external programs
                             // when multiple programs launched by mpirun
                             // set by -mpicolor command line arg
  static bool has_git_info();
  static const char *git_commit();
  static const char *git_branch();
  static const char *git_descriptor();

  LAMMPS(int, char **, MPI_Comm);
  ~LAMMPS();
  void create();
  void post_create();
  void init();
  void destroy();

 private:
  /// Default constructor. Declared private to prohibit its use
  LAMMPS(){};
  /// Copy constructor. Declared private to prohibit its use
  LAMMPS(const LAMMPS &){};
};

}    // namespace LAMMPS_NS

#endif
