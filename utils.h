#ifndef LMP_UTILS_H
#define LMP_UTILS_H
#include "fmt/format.h"
#include "lmptype.h"
#include <mpi.h>
#include <vector>
namespace LAMMPS_NS {
class Error;
class LAMMPS;
namespace utils {
bool strmatch(const std::string &text, const std::string &pattern);
std::string strfind(const std::string &text, const std::string &pattern);
void missing_cmd_args(const std::string &file, int line, const std::string &cmd,
                      Error *error);
void fmtargs_logmesg(LAMMPS *lmp, fmt::string_view format,
                     fmt::format_args args);
template <typename... Args>
void logmesg(LAMMPS *lmp, const std::string &format, Args &&... args) {
  fmtargs_logmesg(lmp, format, fmt::make_format_args(args...));
}
void logmesg(LAMMPS *lmp, const std::string &mesg);
std::string errorurl(int errorcode);
void flush_buffers(LAMMPS *lmp);
std::string getsyserror();
char *fgets_trunc(char *s, int size, FILE *fp);
void sfgets(const char *srcname, int srcline, char *s, int size, FILE *fp,
            const char *filename, Error *error);
int read_lines_from_file(FILE *fp, int nlines, int nmax, char *buffer, int me,
                         MPI_Comm comm);
int logical(const char *file, int line, const std::string &str, bool do_abort,
            LAMMPS *lmp);
int logical(const char *file, int line, const char *str, bool do_abort,
            LAMMPS *lmp);
double numeric(const char *file, int line, const std::string &str,
               bool do_abort, LAMMPS *lmp);
double numeric(const char *file, int line, const char *str, bool do_abort,
               LAMMPS *lmp);
int inumeric(const char *file, int line, const std::string &str, bool do_abort,
             LAMMPS *lmp);
int inumeric(const char *file, int line, const char *str, bool do_abort,
             LAMMPS *lmp);
bigint bnumeric(const char *file, int line, const std::string &str,
                bool do_abort, LAMMPS *lmp);
bigint bnumeric(const char *file, int line, const char *str, bool do_abort,
                LAMMPS *lmp);
tagint tnumeric(const char *file, int line, const std::string &str,
                bool do_abort, LAMMPS *lmp);
tagint tnumeric(const char *file, int line, const char *str, bool do_abort,
                LAMMPS *lmp);
template <typename TYPE>
void bounds(const char *file, int line, const std::string &str, bigint nmin,
            bigint nmax, TYPE &nlo, TYPE &nhi, Error *error);
char *strdup(const std::string &text);
std::string lowercase(const std::string &line);
std::string uppercase(const std::string &line);
std::string trim(const std::string &line);
std::string trim_comment(const std::string &line);
std::string star_subst(const std::string &name, bigint step, int pad);
std::string strip_style_suffix(const std::string &style, LAMMPS *lmp);
inline bool has_utf8(const std::string &line) {
  for (auto c : line)
    if (c & 0x80U)
      return true;
  return false;
}
std::string utf8_subst(const std::string &line);
size_t count_words(const std::string &text, const std::string &separators);
size_t count_words(const std::string &text);
size_t count_words(const char *text);
size_t trim_and_count_words(const std::string &text,
                            const std::string &separators = " \t\r\n\f");
std::string join_words(const std::vector<std::string> &words,
                       const std::string &sep);
std::vector<std::string> split_words(const std::string &text);
bool is_integer(const std::string &str);
bool is_double(const std::string &str);
bool is_id(const std::string &str);
int is_type(const std::string &str);
enum { NOCONVERT = 0, METAL2REAL = 1, REAL2METAL = 1 << 1 };
enum { UNKNOWN = 0, ENERGY };
int get_supported_conversions(const int property);
double get_conversion_factor(const int property, const int conversion);
int binary_search(const double needle, const int n, const double *haystack);
void merge_sort(int *index, int num, void *ptr, int (*comp)(int, int, void *));
} // namespace utils
} // namespace LAMMPS_NS
#endif
