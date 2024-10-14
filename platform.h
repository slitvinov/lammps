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
void usleep(int usec);
std::string os_info();
std::string cxx_standard();
std::string compiler_info();
std::string openmp_standard();
std::string mpi_vendor();
std::string mpi_info(int &major, int &minor);
std::string compress_info();
int putenv(const std::string &vardef);
int unsetenv(const std::string &variable);
std::vector<std::string> list_pathenv(const std::string &var);
void *dlopen(const std::string &fname);
std::string dlerror();
int dlclose(void *handle);
void *dlsym(void *handle, const std::string &symbol);
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
const char *guesspath(FILE *fp, char *buf, int len);
bool is_console(FILE *fp);
std::string current_directory();
bool path_is_directory(const std::string &path);
std::vector<std::string> list_directory(const std::string &dir);
std::string find_exe_path(const std::string &cmd);
int chdir(const std::string &path);
int mkdir(const std::string &path);
int rmdir(const std::string &path);
int unlink(const std::string &path);
bigint ftell(FILE *fp);
constexpr bigint END_OF_FILE = -1;
int fseek(FILE *fp, bigint pos);
int ftruncate(FILE *fp, bigint length);
FILE *popen(const std::string &cmd, const std::string &mode);
int pclose(FILE *fp);
std::string path_basename(const std::string &path);
std::string path_dirname(const std::string &path);
std::string path_join(const std::string &a, const std::string &b);
bool file_is_readable(const std::string &path);
bool has_compress_extension(const std::string &file);
FILE *compressed_read(const std::string &file);
FILE *compressed_write(const std::string &file);
} // namespace platform
} // namespace LAMMPS_NS
#endif
