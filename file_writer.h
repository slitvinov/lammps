#ifndef LMP_FILE_WRITER_H
#define LMP_FILE_WRITER_H
#include <string>
namespace LAMMPS_NS {
class FileWriter {
public:
  FileWriter() = default;
  virtual ~FileWriter() = default;
  virtual void open(const std::string &path, bool append = false) = 0;
  virtual void close() = 0;
  virtual void flush() = 0;
  virtual size_t write(const void *buffer, size_t length) = 0;
  virtual bool isopen() const = 0;
};
class FileWriterException : public std::exception {
  std::string message;

public:
  FileWriterException(const std::string &msg) : message(msg) {}
  const char *what() const noexcept override { return message.c_str(); }
};
} // namespace LAMMPS_NS
#endif
