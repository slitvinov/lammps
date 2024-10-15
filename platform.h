#ifndef LMP_PLATFORM_H
#define LMP_PLATFORM_H
#include "lmptype.h"
#include <cstdio>
#include <string>
#include <vector>
namespace LAMMPS_NS {
namespace platform {
double walltime();
#if !defined(_WIN32)
constexpr char filepathsep[] = "/";
#else
constexpr char filepathsep[] = "\\/";
#endif
#if !defined(_WIN32)
constexpr char pathvarsep = ':';
#else
constexpr char pathvarsep = ';';
#endif
constexpr bigint END_OF_FILE = -1;
} // namespace platform
} // namespace LAMMPS_NS
#endif
