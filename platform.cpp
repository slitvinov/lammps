#include "platform.h"
#include "fmt/format.h"
#include "utils.h"
#include <chrono>
#include <cstring>
#include <deque>
#include <dirent.h>
#include <dlfcn.h>
#include <exception>
#include <mpi.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <thread>
#include <unistd.h>
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
static auto initial_time = std::chrono::steady_clock::now();
