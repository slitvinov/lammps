#include "platform.h"
#include "fmt/format.h"
#include "text_file_reader.h"
#include "utils.h"
#include <deque>
#include <exception>
#include <mpi.h>
#if defined(_WIN32)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN 
#endif
#if defined(_WIN32_WINNT)
#undef _WIN32_WINNT
#endif
#define _WIN32_WINNT _WIN32_WINNT_WIN7
#define PSAPI_VERSION 2
#include <direct.h>
#include <io.h>
#include <sys/stat.h>
#include <windows.h>
#else
#include <dirent.h>
#include <dlfcn.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <unistd.h>
#endif
#if defined(__APPLE__)
#include <fcntl.h>
#include <sys/syslimits.h>
#endif
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
    {"lzma", "xz", " --format=lzma > ", " --format=lzma -cdf ", compress_info::LZMA},
    {"lz4", "lz4", " > ", " -cdf ", compress_info::LZ4},
};
static const compress_info &find_compress_type(const std::string &file)
{
  std::size_t dot = file.find_last_of('.');
  if (dot != std::string::npos) {
    const std::string ext = file.substr(dot + 1);
    for (const auto &i : compress_styles) {
      if (i.extension == ext) return i;
    }
  }
  return compress_styles[0];
}
static auto initial_time = std::chrono::steady_clock::now();
using namespace LAMMPS_NS;
#if defined(__clang__)
[[clang::optnone]]
#elif defined(_MSC_VER)
#pragma optimize("",off)
#endif
double platform::cputime()
{
  double rv = 0.0;
#ifdef _WIN32
  FILETIME ct, et, kt, ut;
  union {
    FILETIME ft;
    uint64_t ui;
  } cpu;
  if (GetProcessTimes(GetCurrentProcess(), &ct, &et, &kt, &ut)) {
    cpu.ft = ut;
    rv = cpu.ui * 0.0000001;
  }
#else
  struct rusage ru;
  if (getrusage(RUSAGE_SELF, &ru) == 0) {
    rv = (double) ru.ru_utime.tv_sec;
    rv += (double) ru.ru_utime.tv_usec * 0.000001;
  }
#endif
  return rv;
}
#if defined(__clang__)
#elif defined(_MSC_VER)
#pragma optimize("", on)
#endif
double platform::walltime()
{
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - initial_time).count();
}
void platform::usleep(int usec)
{
  return std::this_thread::sleep_for(std::chrono::microseconds(usec));
}
std::string platform::os_info()
{
  std::string buf;
#if defined(_WIN32)
  char value[1024];
  DWORD value_length = 1024;
  const char *subkey = "SOFTWARE\\Microsoft\\Windows NT\\CurrentVersion";
  const char *entry = "CurrentBuild";
  RegGetValue(HKEY_LOCAL_MACHINE, subkey, entry, RRF_RT_REG_SZ, nullptr, &value,
              (LPDWORD) &value_length);
  value[1023] = '\0';
  auto build = std::string(value);
  if (build == "6002") {
    buf = "Windows Vista";
  } else if (build == "6003") {
    buf = "Windows Server 2008";
  } else if (build == "7601") {
    buf = "Windows 7";
  } else if (build == "9200") {
    buf = "Windows 8";
  } else if (build == "9600") {
    buf = "Windows 8.1";
  } else if (build == "10240") {
    buf = "Windows 10 1507";
  } else if (build == "10586") {
    buf = "Windows 10 1511";
  } else if (build == "14393") {
    buf = "Windows 10 1607";
  } else if (build == "15063") {
    buf = "Windows 10 1703";
  } else if (build == "16299") {
    buf = "Windows 10 1709";
  } else if (build == "17134") {
    buf = "Windows 10 1803";
  } else if (build == "17763") {
    buf = "Windows 10 1809";
  } else if (build == "18362") {
    buf = "Windows 10 1903";
  } else if (build == "18363") {
    buf = "Windows 10 1909";
  } else if (build == "19041") {
    buf = "Windows 10 2004";
  } else if (build == "19042") {
    buf = "Windows 10 20H2";
  } else if (build == "19043") {
    buf = "Windows 10 21H1";
  } else if (build == "19044") {
    buf = "Windows 10 21H2";
  } else if (build == "19045") {
    buf = "Windows 10 22H2";
  } else if (build == "20348") {
    buf = "Windows Server 2022";
  } else if (build == "22000") {
    buf = "Windows 11 21H2";
  } else if (build == "22621") {
    buf = "Windows 11 22H2";
  } else {
    const char *entry = "ProductName";
    RegGetValue(HKEY_LOCAL_MACHINE, subkey, entry, RRF_RT_REG_SZ, nullptr, &value,
                (LPDWORD) &value_length);
    value[1023] = '\0';
    buf = value;
  }
  DWORD fullversion, majorv, minorv, buildv = 0;
  fullversion = GetVersion();
  majorv = (DWORD) (LOBYTE(LOWORD(fullversion)));
  minorv = (DWORD) (HIBYTE(LOWORD(fullversion)));
  if (fullversion < 0x80000000) buildv = (DWORD) (HIWORD(fullversion));
  buf += ", Windows ABI " + std::to_string(majorv) + "." + std::to_string(minorv) + " (" +
      std::to_string(buildv) + ") on ";
  SYSTEM_INFO si;
  GetSystemInfo(&si);
  switch (si.wProcessorArchitecture) {
    case PROCESSOR_ARCHITECTURE_AMD64:
      buf += "x86_64";
      break;
    case PROCESSOR_ARCHITECTURE_ARM:
      buf += "arm";
      break;
    case PROCESSOR_ARCHITECTURE_IA64:
      buf += "ia64";
      break;
    case PROCESSOR_ARCHITECTURE_INTEL:
      buf += "i386";
      break;
    default:
      buf += "(unknown)";
  }
#else
  struct utsname ut;
  uname(&ut);
  buf = ut.sysname;
  if (platform::file_is_readable("/etc/os-release")) {
    try {
      TextFileReader reader("/etc/os-release", "");
      while (true) {
        auto words = reader.next_values(0, "=");
        if ((words.count() > 1) && (words.next_string() == "PRETTY_NAME")) {
          buf += " " + utils::trim(words.next_string());
          break;
        }
      }
    } catch (std::exception &e) {
      ;
    }
  }
  buf += std::string(" ") + ut.release + " " + ut.machine;
#endif
  return buf;
}
std::string platform::cxx_standard()
{
#if __cplusplus > 202002L
  return "newer than C++20";
#elif __cplusplus == 202002L
  return "C++20";
#elif __cplusplus == 201703L
  return "C++17";
#elif __cplusplus == 201402L
  return "C++14";
#elif __cplusplus == 201103L
  return "C++11";
#elif __cplusplus == 199711L
  return "C++98";
#else
  return "unknown";
#endif
}
std::string platform::compiler_info()
{
  std::string buf = "(Unknown)";
#if defined(__INTEL_LLVM_COMPILER)
  double version = static_cast<double>(__INTEL_LLVM_COMPILER) * 0.01;
  buf = fmt::format("Intel LLVM C++ {:.1f} / {}", version, __VERSION__);
#elif defined(__ibmxl__)
  buf = fmt::format("IBM XL C/C++ (Clang) {}.{}.{}", __ibmxl_version__, __ibmxl_release__,
                    __ibmxl_modification__);
#elif defined(__clang__)
  buf = fmt::format("Clang C++ {}", __VERSION__);
#elif defined(__PGI)
  buf = fmt::format("PGI C++ {}.{}", __PGIC__, __PGIC_MINOR__);
#elif defined(__INTEL_COMPILER)
#if !defined(__VERSION__)
#define __VERSION__ __INTEL_COMPILER_BUILD_DATE
#endif
  double version = static_cast<double>(__INTEL_COMPILER) * 0.01;
  buf = fmt::format("Intel Classic C++ {:.2f}.{} / {}", version, __INTEL_COMPILER_UPDATE,
                    __VERSION__);
#elif defined(__MINGW64__)
  buf = fmt::format("MinGW-w64 64bit {}.{} / GNU C++ {}", __MINGW64_VERSION_MAJOR,
                    __MINGW64_VERSION_MINOR, __VERSION__);
#elif defined(__MINGW32__)
  buf = fmt::format("MinGW-w64 32bit {}.{} / GNU C++ {}", __MINGW32_MAJOR_VERSION,
                    __MINGW32_MINOR_VERSION, __VERSION__);
#elif defined(__GNUC__)
  buf = fmt::format("GNU C++ {}", __VERSION__);
#elif defined(_MSC_VER) && (_MSC_VER >= 1920) && (_MSC_VER < 1930)
  constexpr int major = _MSC_VER / 100;
  constexpr int minor = _MSC_VER - major * 100;
  constexpr int patch = minor - 20;
  buf = fmt::format("Microsoft Visual Studio 2019 Version 16.{}, C/C++ {}.{}", patch, major - 5,
                    minor);
#elif defined(_MSC_VER) && (_MSC_VER >= 1930) && (_MSC_VER < 2000)
  constexpr int major = _MSC_VER / 100;
  constexpr int minor = _MSC_VER - major * 100;
  constexpr int patch = minor - 30;
  buf = fmt::format("Microsoft Visual Studio 2022 Version 17.{}, C/C++ {}.{}", patch, major - 5,
                    minor);
#else
  buf = "(Unknown)";
#endif
  return buf;
}
std::string platform::openmp_standard()
{
#if !defined(_OPENMP)
  return "OpenMP not enabled";
#else
#if _OPENMP > 202011
  return "OpenMP newer than version 5.1";
#elif _OPENMP == 202011
  return "OpenMP 5.1";
#elif _OPENMP == 201811
  return "OpenMP 5.0";
#elif _OPENMP == 201611
  return "OpenMP 5.0 preview 1";
#elif _OPENMP == 201511
  return "OpenMP 4.5";
#elif _OPENMP == 201307
  return "OpenMP 4.0";
#elif _OPENMP == 201107
  return "OpenMP 3.1";
#elif _OPENMP == 200805
  return "OpenMP 3.0";
#elif _OPENMP == 200505
  return "OpenMP 2.5";
#elif _OPENMP == 200203
  return "OpenMP 2.0";
#else
  return "unknown OpenMP version";
#endif
#endif
}
std::string platform::mpi_vendor()
{
#if defined(MPI_STUBS)
  return "MPI STUBS";
#elif defined(OPEN_MPI)
  return "Open MPI";
#elif defined(MPICH_NAME)
  return "MPICH";
#elif defined(I_MPI_VERSION)
  return "Intel MPI";
#elif defined(PLATFORM_MPI)
  return "Platform MPI";
#elif defined(HP_MPI)
  return "HP MPI";
#elif defined(MSMPI_VER)
  char value[1024];
  DWORD value_length = 1024;
  const char *subkey = "SOFTWARE\\Microsoft\\MPI";
  const char *entry = "Version";
  auto rv = RegGetValueA(HKEY_LOCAL_MACHINE, subkey, entry, RRF_RT_REG_SZ, nullptr, &value,
                         (LPDWORD) &value_length);
  std::string buf = "Microsoft MPI";
  if (rv == ERROR_SUCCESS) buf += std::string(" v") + value;
  return buf;
#else
  return "Unknown MPI implementation";
#endif
}
std::string platform::mpi_info(int &major, int &minor)
{
#if (defined(MPI_VERSION) && (MPI_VERSION > 2)) || defined(MPI_STUBS)
  int len = 0;
  static char version[MPI_MAX_LIBRARY_VERSION_STRING];
  MPI_Get_library_version(version, &len);
  if (len > 80) {
    char *ptr = strchr(version + 80, '\n');
    if (ptr) *ptr = '\0';
  }
#else
  constexpr int MAX_VERSION_STRING = 32;
  static char version[MAX_VERSION_STRING];
  strncpy(version, mpi_vendor().c_str(), MAX_VERSION_STRING);
#endif
#if defined(MPI_VERSION)
  MPI_Get_version(&major, &minor);
#else
  major = 1;
  minor = 0;
#endif
  return {version};
}
std::string platform::compress_info()
{
  std::string buf = "Available compression formats:\n\n";
  bool none_found = true;
  for (const auto &cmpi : compress_styles) {
    if (cmpi.style == ::compress_info::NONE) continue;
    if (find_exe_path(cmpi.command).size()) {
      none_found = false;
      buf += fmt::format("Extension: .{:6} Command: {}\n", cmpi.extension, cmpi.command);
    }
  }
  if (none_found) buf += "None\n";
  return buf;
}
int platform::putenv(const std::string &vardef)
{
  if (vardef.size() == 0) return -1;
  auto found = vardef.find_first_of('=');
#ifdef _WIN32
  if (found == std::string::npos)
    return _putenv_s(vardef.c_str(), "1");
  else
    return _putenv_s(vardef.substr(0, found).c_str(), vardef.substr(found + 1).c_str());
#else
  if (found == std::string::npos)
    return setenv(vardef.c_str(), "", 1);
  else
    return setenv(vardef.substr(0, found).c_str(), vardef.substr(found + 1).c_str(), 1);
#endif
  return -1;
}
int platform::unsetenv(const std::string &variable)
{
  if (variable.size() == 0) return -1;
#ifdef _WIN32
  const char *ptr = getenv(variable.c_str());
  if (!ptr) return -1;
  return _putenv_s(variable.c_str(), "");
#else
  return ::unsetenv(variable.c_str());
#endif
}
std::vector<std::string> platform::list_pathenv(const std::string &var)
{
  std::vector<std::string> dirs;
  const char *ptr = getenv(var.c_str());
  if (ptr == nullptr) return dirs;
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
std::string platform::find_exe_path(const std::string &cmd)
{
  if (cmd.size() == 0) return "";
  auto pathdirs = list_pathenv("PATH");
#ifdef _WIN32
  pathdirs.insert(pathdirs.begin(), ".");
#else
  struct stat info;
#endif
  for (const auto &dir : pathdirs) {
    std::string exe = path_join(dir, cmd);
#ifdef _WIN32
    const char *extensions[] = {".exe", ".com", ".bat", nullptr};
    for (auto ext = extensions; *ext != nullptr; ++ext) {
      auto exe_path = exe + *ext;
      if (file_is_readable(exe_path)) return exe_path;
    }
#else
    memset(&info, 0, sizeof(info));
    if (stat(exe.c_str(), &info) != 0) continue;
    if ((info.st_mode & (S_IXOTH | S_IXGRP | S_IXUSR)) != 0) return exe;
#endif
  }
  return "";
}
#ifdef _WIN32
void *platform::dlopen(const std::string &fname)
{
  return (void *) LoadLibrary(fname.c_str());
}
std::string platform::dlerror()
{
  return "";
}
int platform::dlclose(void *handle)
{
  return (FreeLibrary((HINSTANCE) handle) == 0);
}
void *platform::dlsym(void *handle, const std::string &symbol)
{
  return (void *) GetProcAddress((HINSTANCE) handle, symbol.c_str());
}
#else
void *platform::dlopen(const std::string &fname)
{
  return ::dlopen(fname.c_str(), RTLD_NOW | RTLD_GLOBAL);
}
std::string platform::dlerror()
{
  const char *errmesg = ::dlerror();
  if (errmesg)
    return {errmesg};
  else
    return {""};
}
int platform::dlclose(void *handle)
{
  return ::dlclose(handle);
}
void *platform::dlsym(void *handle, const std::string &symbol)
{
  return ::dlsym(handle, symbol.c_str());
}
#endif
const char *platform::guesspath(FILE *fp, char *buf, int len)
{
  if ((buf == nullptr) || (len < 16)) return nullptr;
  memset(buf, 0, len);
  len--;
#if defined(__linux__)
  int fd = fileno(fp);
  if (readlink((std::string("/proc/self/fd/") + std::to_string(fd)).c_str(), buf, len) <= 0)
    strncpy(buf, "(unknown)", len);
#elif defined(__APPLE__)
  int fd = fileno(fp);
  char filepath[PATH_MAX];
  if (fcntl(fd, F_GETPATH, filepath) != -1)
    strncpy(buf, filepath, len);
  else
    strncpy(buf, "(unknown)", len);
#elif defined(_WIN32)
  char filepath[MAX_PATH];
  HANDLE h = (HANDLE) _get_osfhandle(_fileno(fp));
  if (GetFinalPathNameByHandleA(h, filepath, MAX_PATH, FILE_NAME_NORMALIZED) > 0)
    strncpy(buf, filepath, len);
  else
    strncpy(buf, "(unknown)", len);
#else
  strncpy(buf, "(unknown)", len);
#endif
  return buf;
}
bool platform::is_console(FILE *fp)
{
  if (!fp) return false;
#if defined(_WIN32)
  return (_isatty(_fileno(fp)) == 1);
#else
  return (isatty(fileno(fp)) == 1);
#endif
}
std::string platform::current_directory()
{
  std::string cwd;
#if defined(_WIN32)
  char *buf = new char[MAX_PATH];
  if (_getcwd(buf, MAX_PATH)) { cwd = buf; }
  delete[] buf;
#else
  auto buf = new char[PATH_MAX];
  if (::getcwd(buf, PATH_MAX)) { cwd = buf; }
  delete[] buf;
#endif
  return cwd;
}
bool platform::path_is_directory(const std::string &path)
{
#if defined(_WIN32)
  struct _stat info;
  memset(&info, 0, sizeof(info));
  if (_stat(path.c_str(), &info) != 0) return false;
#else
  struct stat info;
  memset(&info, 0, sizeof(info));
  if (stat(path.c_str(), &info) != 0) return false;
#endif
  return ((info.st_mode & S_IFDIR) != 0);
}
std::vector<std::string> platform::list_directory(const std::string &dir)
{
  std::vector<std::string> files;
  if (!path_is_directory(dir)) return files;
#if defined(_WIN32)
  HANDLE handle;
  WIN32_FIND_DATA fd;
  std::string searchname = dir + filepathsep[0] + "*";
  handle = FindFirstFile(searchname.c_str(), &fd);
  if (handle == ((HANDLE) -1)) return files;
  while (FindNextFile(handle, &fd)) {
    std::string entry(fd.cFileName);
    if ((entry == "..") || (entry == ".")) continue;
    files.push_back(entry);
  }
  FindClose(handle);
#else
  std::string dirname = dir + filepathsep[0];
  DIR *handle = opendir(dirname.c_str());
  if (handle == nullptr) return files;
  struct dirent *fd;
  while ((fd = readdir(handle)) != nullptr) {
    std::string entry(fd->d_name);
    if ((entry == "..") || (entry == ".")) continue;
    files.push_back(entry);
  }
  closedir(handle);
#endif
  return files;
}
int platform::chdir(const std::string &path)
{
#if defined(_WIN32)
  return ::_chdir(path.c_str());
#else
  return ::chdir(path.c_str());
#endif
}
int platform::mkdir(const std::string &path)
{
  std::deque<std::string> dirlist = {path};
  std::string dirname = path_dirname(path);
  while ((dirname != ".") && (dirname != "")) {
    dirlist.push_front(dirname);
    dirname = path_dirname(dirname);
  }
  int rv;
  for (const auto &dir : dirlist) {
    if (!path_is_directory(dir)) {
#if defined(_WIN32)
      rv = ::_mkdir(dir.c_str());
#else
      rv = ::mkdir(dir.c_str(), S_IRWXU | S_IRGRP | S_IXGRP);
#endif
      if (rv != 0) return rv;
    }
  }
  return 0;
}
int platform::rmdir(const std::string &path)
{
  auto entries = list_directory(path);
  for (const auto &entry : entries) {
    const auto newpath = path_join(path, entry);
    if (path_is_directory(newpath))
      rmdir(newpath);
    else
      unlink(newpath);
  }
#if defined(_WIN32)
  return ::_rmdir(path.c_str());
#else
  return ::rmdir(path.c_str());
#endif
}
int platform::unlink(const std::string &path)
{
#if defined(_WIN32)
  return ::_unlink(path.c_str());
#else
  return ::unlink(path.c_str());
#endif
}
bigint platform::ftell(FILE *fp)
{
#if defined(_WIN32)
  return (bigint)::_ftelli64(fp);
#else
  return (bigint)::ftell(fp);
#endif
}
int platform::fseek(FILE *fp, bigint pos)
{
#if defined(_WIN32)
  if (pos == platform::END_OF_FILE)
    return ::_fseeki64(fp, 0, SEEK_END);
  else
    return ::_fseeki64(fp, (__int64) pos, SEEK_SET);
#else
  if (pos == platform::END_OF_FILE)
    return ::fseek(fp, 0, SEEK_END);
  else
    return ::fseek(fp, (long) pos, SEEK_SET);
#endif
}
int platform::ftruncate(FILE *fp, bigint length)
{
#if defined(_WIN32)
  HANDLE h = (HANDLE) _get_osfhandle(_fileno(fp));
  LARGE_INTEGER li_start, li_length;
  li_start.QuadPart = (int64_t) 0;
  li_length.QuadPart = (int64_t) length;
  if (SetFilePointerEx(h, li_start, NULL, FILE_CURRENT) &&
      SetFilePointerEx(h, li_length, NULL, FILE_BEGIN) && SetEndOfFile(h)) {
    return 0;
  } else {
    return 1;
  }
#else
  platform::fseek(fp, length);
  return ::ftruncate(fileno(fp), (off_t) length);
#endif
}
FILE *platform::popen(const std::string &cmd, const std::string &mode)
{
  FILE *fp = nullptr;
#if defined(_WIN32)
  if (mode == "r")
    fp = ::_popen(cmd.c_str(), "rb");
  else if (mode == "w")
    fp = ::_popen(cmd.c_str(), "wb");
#else
  if (mode == "r")
    fp = ::popen(cmd.c_str(), "r");
  else if (mode == "w")
    fp = ::popen(cmd.c_str(), "w");
#endif
  return fp;
}
int platform::pclose(FILE *fp)
{
#if defined(_WIN32)
  return ::_pclose(fp);
#else
  return ::pclose(fp);
#endif
}
std::string platform::path_basename(const std::string &path)
{
  size_t start = path.find_last_of(platform::filepathsep);
  if (start == std::string::npos) {
    start = 0;
  } else {
    start += 1;
  }
  return path.substr(start);
}
std::string platform::path_dirname(const std::string &path)
{
  size_t start = path.find_last_of(platform::filepathsep);
  if (start == std::string::npos) return ".";
  return path.substr(0, start);
}
std::string platform::path_join(const std::string &a, const std::string &b)
{
  if (a.empty()) return b;
  if (b.empty()) return a;
  std::string joined = a;
  while (joined.find_last_of(platform::filepathsep) == joined.size() - 1) {
    for (const auto &s : platform::filepathsep)
      if (joined.back() == s) joined.pop_back();
  }
  std::size_t skip = 0;
  while (b.find_first_of(platform::filepathsep, skip) == skip) ++skip;
  joined += platform::filepathsep[0] + b.substr(skip);
  return joined;
}
bool platform::file_is_readable(const std::string &path)
{
  FILE *fp = fopen(path.c_str(), "r");
  if (fp) {
    fclose(fp);
    return true;
  }
  return false;
}
bool platform::has_compress_extension(const std::string &file)
{
  return find_compress_type(file).style != ::compress_info::NONE;
}
FILE *platform::compressed_read(const std::string &file)
{
  FILE *fp = nullptr;
#if defined(LAMMPS_GZIP)
  auto compress = find_compress_type(file);
  if (compress.style == ::compress_info::NONE) return nullptr;
  if (find_exe_path(compress.command).size())
    fp = popen((compress.command + compress.uncompressflags + "\"" + file + "\""), "r");
#endif
  return fp;
}
FILE *platform::compressed_write(const std::string &file)
{
  FILE *fp = nullptr;
#if defined(LAMMPS_GZIP)
  auto compress = find_compress_type(file);
  if (compress.style == ::compress_info::NONE) return nullptr;
  if (find_exe_path(compress.command).size())
    fp = popen((compress.command + compress.compressflags + "\"" + file + "\""), "w");
#endif
  return fp;
}
