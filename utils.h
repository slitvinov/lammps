#ifndef LMP_UTILS_H
#define LMP_UTILS_H
#include "lmptype.h"
#include <mpi.h>
#include <vector>
namespace LAMMPS_NS {
class Error;
class LAMMPS;
namespace utils {
bool strmatch(const std::string &text, const std::string &pattern);
void sfgets(const char *srcname, int srcline, char *s, int size, FILE *fp,
            const char *filename, Error *error);
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
template <typename TYPE>
void bounds(const char *file, int line, const std::string &str, bigint nmin,
            bigint nmax, TYPE &nlo, TYPE &nhi, Error *error);
char *strdup(const std::string &text);
std::string trim(const std::string &line);
std::string strip_style_suffix(const std::string &style, LAMMPS *lmp);
std::vector<std::string> split_words(const std::string &text);
bool is_integer(const std::string &str);
bool is_double(const std::string &str);
bool is_id(const std::string &str);
enum { NOCONVERT = 0, METAL2REAL = 1, REAL2METAL = 1 << 1 };
enum { UNKNOWN = 0, ENERGY };
int get_supported_conversions(const int property);
double get_conversion_factor(const int property, const int conversion);
} // namespace utils
} // namespace LAMMPS_NS
#endif
