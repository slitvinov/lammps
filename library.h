/* -*- c -*- ------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LAMMPS_LIBRARY_H
#define LAMMPS_LIBRARY_H

/* C style library interface to LAMMPS which allows to create and
 * control instances of the LAMMPS C++ class and exchange data with it.
 * The C bindings are the basis for the Python and Fortran modules.
 *
 * If needed, new LAMMPS-specific functions can be added to expose
 * additional LAMMPS functionality to this library interface. */

/* We follow the behavior of regular LAMMPS compilation and assume
 * -DLAMMPS_SMALLBIG when no define is set. */

#if !defined(LAMMPS_BIGBIG) && !defined(LAMMPS_SMALLBIG) && !defined(LAMMPS_SMALLSMALL)
#define LAMMPS_SMALLBIG
#endif

/* To allow including the library interface without MPI */

#if defined(LAMMPS_LIB_MPI)
#include <mpi.h>
#endif

#if defined(LAMMPS_BIGBIG) || defined(LAMMPS_SMALLBIG)
#include <stdint.h> /* for int64_t */
#endif

/** Data type constants for extracting data from atoms, computes and fixes
 *
 * Must be kept in sync with the equivalent constants in ``python/lammps/constants.py``,
 * ``fortran/lammps.f90``, ``tools/swig/lammps.i``, and
 *``examples/COUPLE/plugin/liblammpsplugin.h`` */

enum _LMP_DATATYPE_CONST {
  LAMMPS_INT = 0,       /*!< 32-bit integer (array) */
  LAMMPS_INT_2D = 1,    /*!< two-dimensional 32-bit integer array */
  LAMMPS_DOUBLE = 2,    /*!< 64-bit double (array) */
  LAMMPS_DOUBLE_2D = 3, /*!< two-dimensional 64-bit double array */
  LAMMPS_INT64 = 4,     /*!< 64-bit integer (array) */
  LAMMPS_INT64_2D = 5,  /*!< two-dimensional 64-bit integer array */
  LAMMPS_STRING = 6     /*!< C-String */
};

/** Style constants for extracting data from computes and fixes.
 *
 * Must be kept in sync with the equivalent constants in ``python/lammps/constants.py``,
 * ``fortran/lammps.f90``, ``tools/swig/lammps.i``, and
 * ``examples/COUPLE/plugin/liblammpsplugin.h`` */

enum _LMP_STYLE_CONST {
  LMP_STYLE_GLOBAL = 0, /*!< return global data */
  LMP_STYLE_ATOM = 1,   /*!< return per-atom data */
  LMP_STYLE_LOCAL = 2   /*!< return local data */
};

/** Type and size constants for extracting data from computes and fixes.
 *
 * Must be kept in sync with the equivalent constants in ``python/lammps/constants.py``,
 * ``fortran/lammps.f90``, ``tools/swig/lammps.i``, and
 * ``examples/COUPLE/plugin/liblammpsplugin.h`` */

enum _LMP_TYPE_CONST {
  LMP_TYPE_SCALAR = 0, /*!< return scalar */
  LMP_TYPE_VECTOR = 1, /*!< return vector */
  LMP_TYPE_ARRAY = 2,  /*!< return array */
  LMP_SIZE_VECTOR = 3, /*!< return length of vector */
  LMP_SIZE_ROWS = 4,   /*!< return number of rows */
  LMP_SIZE_COLS = 5    /*!< return number of columns */
};

/** Error codes to select the suitable function in the Error class
 *
 * Must be kept in sync with the equivalent constants in ``python/lammps/constants.py``,
 * ``fortran/lammps.f90``, ``tools/swig/lammps.i``, and
 * ``examples/COUPLE/plugin/liblammpsplugin.h`` */

enum _LMP_ERROR_CONST {
  LMP_ERROR_WARNING = 0, /*!< call Error::warning() */
  LMP_ERROR_ONE = 1,     /*!< called from one MPI rank */
  LMP_ERROR_ALL = 2,     /*!< called from all MPI ranks */
  LMP_ERROR_WORLD = 4,   /*!< error on Comm::world */
  LMP_ERROR_UNIVERSE = 8 /*!< error on Comm::universe */
};

/** Variable style constants for extracting data from variables.
 *
 * Must be kept in sync with the equivalent constants in ``python/lammps/constants.py``,
 * ``fortran/lammps.f90``, ``tools/swig/lammps.i``, and
 * ``examples/COUPLE/plugin/liblammpsplugin.h`` */

enum _LMP_VAR_CONST {
  LMP_VAR_EQUAL = 0,  /*!< compatible with equal-style variables */
  LMP_VAR_ATOM = 1,   /*!< compatible with atom-style variables */
  LMP_VAR_VECTOR = 2, /*!< compatible with vector-style variables */
  LMP_VAR_STRING = 3  /*!< return value will be a string (catch-all) */
};

/* Ifdefs to allow this file to be included in C and C++ programs */

#ifdef __cplusplus
extern "C" {
#endif

#endif /* LAMMPS_LIBRARY_H */
