#if !defined(_CRT_SECURE_NO_WARNINGS) && defined(_MSC_VER)
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "fmt/os.h"
#include <climits>
#if FMT_USE_FCNTL
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#endif
namespace {
#if   FMT_USE_FCNTL
using rwresult = ssize_t;
inline std::size_t convert_rwcount(std::size_t count) { return count; }
#endif
} // namespace
FMT_BEGIN_NAMESPACE
buffered_file::~buffered_file() noexcept {
  if (file_ && FMT_SYSTEM(fclose(file_)) != 0)
    report_system_error(errno, "cannot close file");
}
buffered_file::buffered_file(cstring_view filename, cstring_view mode) {
  FMT_RETRY_VAL(file_, FMT_SYSTEM(fopen(filename.c_str(), mode.c_str())),
                nullptr);
  if (!file_)
    FMT_THROW(system_error(errno, "cannot open file {}", filename.c_str()));
}
void buffered_file::close() {
  if (!file_)
    return;
  int result = FMT_SYSTEM(fclose(file_));
  file_ = nullptr;
  if (result != 0)
    FMT_THROW(system_error(errno, "cannot close file"));
}
int buffered_file::descriptor() const {
  int fd = FMT_POSIX_CALL(fileno(file_));
  if (fd == -1)
    FMT_THROW(system_error(errno, "cannot get file descriptor"));
  return fd;
}
#if FMT_USE_FCNTL
file::file(cstring_view path, int oflag) {
  constexpr mode_t mode =
      S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH;
  FMT_RETRY(fd_, FMT_POSIX_CALL(open(path.c_str(), oflag, mode)));
  if (fd_ == -1)
    FMT_THROW(system_error(errno, "cannot open file {}", path.c_str()));
}
file::~file() noexcept {
  if (fd_ != -1 && FMT_POSIX_CALL(close(fd_)) != 0)
    report_system_error(errno, "cannot close file");
}
void file::close() {
  if (fd_ == -1)
    return;
  int result = FMT_POSIX_CALL(close(fd_));
  fd_ = -1;
  if (result != 0)
    FMT_THROW(system_error(errno, "cannot close file"));
}
long long file::size() const {
  using Stat = struct stat;
  Stat file_stat = Stat();
  if (FMT_POSIX_CALL(fstat(fd_, &file_stat)) == -1)
    FMT_THROW(system_error(errno, "cannot get file attributes"));
  static_assert(sizeof(long long) >= sizeof(file_stat.st_size),
                "return type of file::size is not large enough");
  return file_stat.st_size;
}
std::size_t file::read(void *buffer, std::size_t count) {
  rwresult result = 0;
  FMT_RETRY(result, FMT_POSIX_CALL(read(fd_, buffer, convert_rwcount(count))));
  if (result < 0)
    FMT_THROW(system_error(errno, "cannot read from file"));
  return detail::to_unsigned(result);
}
std::size_t file::write(const void *buffer, std::size_t count) {
  rwresult result = 0;
  FMT_RETRY(result, FMT_POSIX_CALL(write(fd_, buffer, convert_rwcount(count))));
  if (result < 0)
    FMT_THROW(system_error(errno, "cannot write to file"));
  return detail::to_unsigned(result);
}
file file::dup(int fd) {
  int new_fd = FMT_POSIX_CALL(dup(fd));
  if (new_fd == -1)
    FMT_THROW(system_error(errno, "cannot duplicate file descriptor {}", fd));
  return file(new_fd);
}
void file::dup2(int fd) {
  int result = 0;
  FMT_RETRY(result, FMT_POSIX_CALL(dup2(fd_, fd)));
  if (result == -1) {
    FMT_THROW(system_error(errno, "cannot duplicate file descriptor {} to {}",
                           fd_, fd));
  }
}
void file::dup2(int fd, std::error_code &ec) noexcept {
  int result = 0;
  FMT_RETRY(result, FMT_POSIX_CALL(dup2(fd_, fd)));
  if (result == -1)
    ec = std::error_code(errno, std::generic_category());
}
void file::pipe(file &read_end, file &write_end) {
  read_end.close();
  write_end.close();
  int fds[2] = {};
  int result = FMT_POSIX_CALL(pipe(fds));
  if (result != 0)
    FMT_THROW(system_error(errno, "cannot create pipe"));
  read_end = file(fds[0]);
  write_end = file(fds[1]);
}
buffered_file file::fdopen(const char *mode) {
#if defined(__MINGW32__) && defined(_POSIX_)
  FILE *f = ::fdopen(fd_, mode);
#else
  FILE *f = FMT_POSIX_CALL(fdopen(fd_, mode));
#endif
  if (!f)
    FMT_THROW(
        system_error(errno, "cannot associate stream with file descriptor"));
  buffered_file bf(f);
  fd_ = -1;
  return bf;
}
long getpagesize() {
  long size = FMT_POSIX_CALL(sysconf(_SC_PAGESIZE));
  if (size < 0)
    FMT_THROW(system_error(errno, "cannot get memory page size"));
  return size;
}
FMT_API void ostream::grow(size_t) {
  if (this->size() == this->capacity())
    flush();
}
#endif
FMT_END_NAMESPACE
