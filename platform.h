#ifndef LMP_PLATFORM_H
#define LMP_PLATFORM_H
#include "lmptype.h"
#include <cstdio>
#include <string>
#include <vector>
namespace LAMMPS_NS {
namespace platform {
double cputime();
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
bool is_console(FILE *fp);
std::string current_directory();
bool path_is_directory(const std::string &path);
constexpr bigint END_OF_FILE = -1;
} // namespace platform
} // namespace LAMMPS_NS
#endif
