#include "text_file_reader.h"
#include "fmt/format.h"
#include "tokenizer.h"
#include "utils.h"
#include <cstring>
#include <utility>
using namespace LAMMPS_NS;
TextFileReader::TextFileReader(const std::string &filename, const std::string &filetype) :
    filetype(filetype), closefp(true), line(nullptr), ignore_comments(true)
{
  set_bufsize(1024);
  fp = fopen(filename.c_str(), "r");
  if (fp == nullptr) {
    throw FileReaderException(
        fmt::format("cannot open {} file {}: {}", filetype, filename, utils::getsyserror()));
  }
}
TextFileReader::TextFileReader(FILE *fp, std::string filetype) :
    filetype(std::move(filetype)), closefp(false), line(nullptr), fp(fp), ignore_comments(true)
{
  set_bufsize(1024);
  if (fp == nullptr) throw FileReaderException("Invalid file descriptor");
}
TextFileReader::~TextFileReader()
{
  if (closefp) fclose(fp);
  delete[] line;
}
void TextFileReader::set_bufsize(int newsize)
{
  if (newsize < 100) {
    throw FileReaderException(
        fmt::format("line buffer size {} for {} file too small, must be > 100", newsize, filetype));
  }
  delete[] line;
  bufsize = newsize;
  line = new char[bufsize];
}
void TextFileReader::rewind()
{
  ::rewind(fp);
}
void TextFileReader::skip_line()
{
  char *ptr = fgets(line, bufsize, fp);
  if (ptr == nullptr) {
    throw EOFException(fmt::format("Missing line in {} file!", filetype));
  }
}
char *TextFileReader::next_line(int nparams)
{
  int n = 0;
  int nwords = 0;
  char *ptr = fgets(line, bufsize, fp);
  if (ptr == nullptr) {
    return nullptr;
  }
  if (ignore_comments && (ptr = strchr(line, '#'))) *ptr = '\0';
  nwords = utils::count_words(line);
  if (nwords > 0) n = strlen(line);
  while (nwords == 0 || nwords < nparams) {
    ptr = fgets(&line[n], bufsize - n, fp);
    if (ptr == nullptr) {
      if (nwords > 0 && nwords < nparams) {
        throw EOFException(fmt::format("Incorrect format in {} file! {}/{} parameters", filetype,
                                       nwords, nparams));
      }
      return nullptr;
    }
    if (ignore_comments && (ptr = strchr(line, '#'))) *ptr = '\0';
    nwords += utils::count_words(&line[n]);
    if (nwords > 0) { n = strlen(line); }
  }
  return line;
}
void TextFileReader::next_dvector(double *list, int n)
{
  int i = 0;
  while (i < n) {
    char *ptr = next_line();
    if (ptr == nullptr) {
      if (i < n) {
        throw FileReaderException(
            fmt::format("Incorrect format in {} file! {}/{} values", filetype, i, n));
      }
    }
    ValueTokenizer values(line);
    while (values.has_next() && i < n) { list[i++] = values.next_double(); }
  }
}
ValueTokenizer TextFileReader::next_values(int nparams, const std::string &separators)
{
  char *ptr = next_line(nparams);
  if (ptr == nullptr) throw EOFException(fmt::format("Missing line in {} file!", filetype));
  return {line, separators};
}
