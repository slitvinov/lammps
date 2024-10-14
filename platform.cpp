#include "platform.h"
#include "fmt/format.h"
#include "utils.h"
#include <deque>
#include <exception>
#include <mpi.h>
#include <dirent.h>
#include <dlfcn.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <unistd.h>
#include <chrono>
#include <cstring>
#include <thread>
struct compress_info {
  enum styles { NONE, GZIP, BZIP2, ZSTD, XZ, LZMA, LZ4 };
  const std::string extension;
  const std::string command;
  const std::string compressflags;
  const std::string uncompressflags;
  const int style;
};
static const std::vector<compress_info> compress_styles = {
    {"", "", "", "", compress_info::NONE},
    {"gz", "gzip", " > ", " -cdf ", compress_info::GZIP},
    {"bz2", "bzip2", " > ", " -cdf ", compress_info::BZIP2},
    {"zst", "zstd", " -q > ", " -cdf ", compress_info::ZSTD},
    {"xz", "xz", " > ", " -cdf ", compress_info::XZ},
    {"lzma", "xz", " --format=lzma > ", " --format=lzma -cdf ",
     compress_info::LZMA},
    {"lz4", "lz4", " > ", " -cdf ", compress_info::LZ4},
};
static const compress_info &find_compress_type(const std::string &file) {
  std::size_t dot = file.find_last_of('.');
  if (dot != std::string::npos) {
    const std::string ext = file.substr(dot + 1);
    for (const auto &i : compress_styles) {
      if (i.extension == ext)
        return i;
    }
  }
  return compress_styles[0];
}
static auto initial_time = std::chrono::steady_clock::now();
using namespace LAMMPS_NS;
#if defined(__clang__)
[[clang::optnone]]
#elif defined(_MSC_VER)
#pragma optimize("", off)
#endif
double
platform::cputime() {
  double rv = 0.0;
  struct rusage ru;
  if (getrusage(RUSAGE_SELF, &ru) == 0) {
    rv = (double)ru.ru_utime.tv_sec;
    rv += (double)ru.ru_utime.tv_usec * 0.000001;
  }
  return rv;
}
#if defined(__clang__)
#elif defined(_MSC_VER)
#pragma optimize("", on)
#endif
double platform::walltime() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                       initial_time)
      .count();
}
void platform::usleep(int usec) {
  return std::this_thread::sleep_for(std::chrono::microseconds(usec));
}
int platform::putenv(const std::string &vardef) {
  if (vardef.size() == 0)
    return -1;
  auto found = vardef.find_first_of('=');
  if (found == std::string::npos)
    return setenv(vardef.c_str(), "", 1);
  else
    return setenv(vardef.substr(0, found).c_str(),
                  vardef.substr(found + 1).c_str(), 1);
  return -1;
}
int platform::unsetenv(const std::string &variable) {
  if (variable.size() == 0)
    return -1;
  return ::unsetenv(variable.c_str());
}
std::vector<std::string> platform::list_pathenv(const std::string &var) {
  std::vector<std::string> dirs;
  const char *ptr = getenv(var.c_str());
  if (ptr == nullptr)
    return dirs;
  std::string pathvar = ptr;
  std::size_t first = 0, next;
  while (true) {
    next = pathvar.find_first_of(pathvarsep, first);
    if (next == std::string::npos) {
      dirs.push_back(pathvar.substr(first));
      break;
    } else {
      dirs.push_back(pathvar.substr(first, next - first));
      first = next + 1;
    }
  }
  return dirs;
}
std::string platform::find_exe_path(const std::string &cmd) {
  if (cmd.size() == 0)
    return "";
  auto pathdirs = list_pathenv("PATH");
  struct stat info;
  for (const auto &dir : pathdirs) {
    std::string exe = path_join(dir, cmd);
    memset(&info, 0, sizeof(info));
    if (stat(exe.c_str(), &info) != 0)
      continue;
    if ((info.st_mode & (S_IXOTH | S_IXGRP | S_IXUSR)) != 0)
      return exe;
  }
  return "";
}
void *platform::dlopen(const std::string &fname) {
  return ::dlopen(fname.c_str(), RTLD_NOW | RTLD_GLOBAL);
}
std::string platform::dlerror() {
  const char *errmesg = ::dlerror();
  if (errmesg)
    return {errmesg};
  else
    return {""};
}
int platform::dlclose(void *handle) { return ::dlclose(handle); }
void *platform::dlsym(void *handle, const std::string &symbol) {
  return ::dlsym(handle, symbol.c_str());
}
const char *platform::guesspath(FILE *fp, char *buf, int len) {
  if ((buf == nullptr) || (len < 16))
    return nullptr;
  memset(buf, 0, len);
  len--;
#if defined(__linux__)
  int fd = fileno(fp);
  if (readlink((std::string("/proc/self/fd/") + std::to_string(fd)).c_str(),
               buf, len) <= 0)
    strncpy(buf, "(unknown)", len);
#else
  strncpy(buf, "(unknown)", len);
#endif
  return buf;
}
bool platform::is_console(FILE *fp) {
  if (!fp)
    return false;
  return (isatty(fileno(fp)) == 1);
}
std::string platform::current_directory() {
  std::string cwd;
  auto buf = new char[PATH_MAX];
  if (::getcwd(buf, PATH_MAX)) {
    cwd = buf;
  }
  delete[] buf;
  return cwd;
}
bool platform::path_is_directory(const std::string &path) {
  struct stat info;
  memset(&info, 0, sizeof(info));
  if (stat(path.c_str(), &info) != 0)
    return false;
  return ((info.st_mode & S_IFDIR) != 0);
}
std::vector<std::string> platform::list_directory(const std::string &dir) {
  std::vector<std::string> files;
  if (!path_is_directory(dir))
    return files;
  std::string dirname = dir + filepathsep[0];
  DIR *handle = opendir(dirname.c_str());
  if (handle == nullptr)
    return files;
  struct dirent *fd;
  while ((fd = readdir(handle)) != nullptr) {
    std::string entry(fd->d_name);
    if ((entry == "..") || (entry == "."))
      continue;
    files.push_back(entry);
  }
  closedir(handle);
  return files;
}
int platform::chdir(const std::string &path) {
  return ::chdir(path.c_str());
}
int platform::mkdir(const std::string &path) {
  std::deque<std::string> dirlist = {path};
  std::string dirname = path_dirname(path);
  while ((dirname != ".") && (dirname != "")) {
    dirlist.push_front(dirname);
    dirname = path_dirname(dirname);
  }
  int rv;
  for (const auto &dir : dirlist) {
    if (!path_is_directory(dir)) {
      rv = ::mkdir(dir.c_str(), S_IRWXU | S_IRGRP | S_IXGRP);
      if (rv != 0)
        return rv;
    }
  }
  return 0;
}
int platform::rmdir(const std::string &path) {
  auto entries = list_directory(path);
  for (const auto &entry : entries) {
    const auto newpath = path_join(path, entry);
    if (path_is_directory(newpath))
      rmdir(newpath);
    else
      unlink(newpath);
  }
  return ::rmdir(path.c_str());
}
int platform::unlink(const std::string &path) {
  return ::unlink(path.c_str());
}
bigint platform::ftell(FILE *fp) {
  return (bigint)::ftell(fp);
}
int platform::fseek(FILE *fp, bigint pos) {
  if (pos == platform::END_OF_FILE)
    return ::fseek(fp, 0, SEEK_END);
  else
    return ::fseek(fp, (long)pos, SEEK_SET);
}
int platform::ftruncate(FILE *fp, bigint length) {
  platform::fseek(fp, length);
  return ::ftruncate(fileno(fp), (off_t)length);
}
FILE *platform::popen(const std::string &cmd, const std::string &mode) {
  FILE *fp = nullptr;
  if (mode == "r")
    fp = ::popen(cmd.c_str(), "r");
  else if (mode == "w")
    fp = ::popen(cmd.c_str(), "w");
  return fp;
}
int platform::pclose(FILE *fp) {
  return ::pclose(fp);
}
std::string platform::path_basename(const std::string &path) {
  size_t start = path.find_last_of(platform::filepathsep);
  if (start == std::string::npos) {
    start = 0;
  } else {
    start += 1;
  }
  return path.substr(start);
}
std::string platform::path_dirname(const std::string &path) {
  size_t start = path.find_last_of(platform::filepathsep);
  if (start == std::string::npos)
    return ".";
  return path.substr(0, start);
}
std::string platform::path_join(const std::string &a, const std::string &b) {
  if (a.empty())
    return b;
  if (b.empty())
    return a;
  std::string joined = a;
  while (joined.find_last_of(platform::filepathsep) == joined.size() - 1) {
    for (const auto &s : platform::filepathsep)
      if (joined.back() == s)
        joined.pop_back();
  }
  std::size_t skip = 0;
  while (b.find_first_of(platform::filepathsep, skip) == skip)
    ++skip;
  joined += platform::filepathsep[0] + b.substr(skip);
  return joined;
}
bool platform::file_is_readable(const std::string &path) {
  FILE *fp = fopen(path.c_str(), "r");
  if (fp) {
    fclose(fp);
    return true;
  }
  return false;
}
bool platform::has_compress_extension(const std::string &file) {
  return find_compress_type(file).style != ::compress_info::NONE;
}
FILE *platform::compressed_read(const std::string &file) {
  FILE *fp = nullptr;
#if defined(LAMMPS_GZIP)
  auto compress = find_compress_type(file);
  if (compress.style == ::compress_info::NONE)
    return nullptr;
  if (find_exe_path(compress.command).size())
    fp = popen(
        (compress.command + compress.uncompressflags + "\"" + file + "\""),
        "r");
#endif
  return fp;
}
FILE *platform::compressed_write(const std::string &file) {
  FILE *fp = nullptr;
#if defined(LAMMPS_GZIP)
  auto compress = find_compress_type(file);
  if (compress.style == ::compress_info::NONE)
    return nullptr;
  if (find_exe_path(compress.command).size())
    fp = popen((compress.command + compress.compressflags + "\"" + file + "\""),
               "w");
#endif
  return fp;
}
