#ifndef LMP_TEXT_FILE_READER_H
#define LMP_TEXT_FILE_READER_H
#include "tokenizer.h"
#include <cstdio>
#include <exception>
#include <string>
namespace LAMMPS_NS {
class TextFileReader {
  std::string filetype;
  bool closefp;
  int bufsize;
  char *line;
  FILE *fp;

public:
  bool ignore_comments;
  TextFileReader(const std::string &filename, const std::string &filetype);
  TextFileReader(FILE *fp, std::string filetype);
  TextFileReader() = delete;
  TextFileReader(const TextFileReader &) = delete;
  TextFileReader(const TextFileReader &&) = delete;
  TextFileReader &operator=(const TextFileReader &) = delete;
  virtual ~TextFileReader();
  void set_bufsize(int);
  void rewind();
  void skip_line();
  char *next_line(int nparams = 0);
  void next_dvector(double *list, int n);
  ValueTokenizer
  next_values(int nparams,
              const std::string &separators = TOKENIZER_DEFAULT_SEPARATORS);
};
class FileReaderException : public std::exception {
  std::string message;

public:
  FileReaderException(const std::string &msg) : message(msg) {}
  const char *what() const noexcept override { return message.c_str(); }
};
class EOFException : public FileReaderException {
public:
  EOFException(const std::string &msg) : FileReaderException(msg) {}
};
} // namespace LAMMPS_NS
#endif
