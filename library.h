#ifndef LAMMPS_LIBRARY_H
#define LAMMPS_LIBRARY_H 
#if !defined(LAMMPS_BIGBIG) && !defined(LAMMPS_SMALLBIG) && !defined(LAMMPS_SMALLSMALL)
#define LAMMPS_SMALLBIG 
#endif
#if defined(LAMMPS_LIB_MPI)
#include <mpi.h>
#endif
#if defined(LAMMPS_BIGBIG) || defined(LAMMPS_SMALLBIG)
#include <stdint.h>
#endif
enum _LMP_DATATYPE_CONST {
  LAMMPS_INT = 0,
  LAMMPS_INT_2D = 1,
  LAMMPS_DOUBLE = 2,
  LAMMPS_DOUBLE_2D = 3,
  LAMMPS_INT64 = 4,
  LAMMPS_INT64_2D = 5,
  LAMMPS_STRING = 6
};
enum _LMP_STYLE_CONST {
  LMP_STYLE_GLOBAL = 0,
  LMP_STYLE_ATOM = 1,
  LMP_STYLE_LOCAL = 2
};
enum _LMP_TYPE_CONST {
  LMP_TYPE_SCALAR = 0,
  LMP_TYPE_VECTOR = 1,
  LMP_TYPE_ARRAY = 2,
  LMP_SIZE_VECTOR = 3,
  LMP_SIZE_ROWS = 4,
  LMP_SIZE_COLS = 5
};
enum _LMP_ERROR_CONST {
  LMP_ERROR_WARNING = 0,
  LMP_ERROR_ONE = 1,
  LMP_ERROR_ALL = 2,
  LMP_ERROR_WORLD = 4,
  LMP_ERROR_UNIVERSE = 8
};
enum _LMP_VAR_CONST {
  LMP_VAR_EQUAL = 0,
  LMP_VAR_ATOM = 1,
  LMP_VAR_VECTOR = 2,
  LMP_VAR_STRING = 3
};
#endif
