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

#include "lammps.h"

#include "input.h"
#include "library.h"


#include <cstdlib>
#include <mpi.h>

#include <stdlib.h>
#include <unistd.h>

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   main program to drive LAMMPS
------------------------------------------------------------------------- */

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  MPI_Comm lammps_comm = MPI_COMM_WORLD;


  try {
    auto lammps = new LAMMPS(argc, argv, lammps_comm);
    lammps->input->file();
    delete lammps;
  } catch (fmt::format_error &fe) {
    fprintf(stderr, "fmt::format_error: %s\n", fe.what());
    MPI_Abort(MPI_COMM_WORLD, 1);
    exit(1);
  }
  MPI_Barrier(lammps_comm);
  MPI_Finalize();
}
