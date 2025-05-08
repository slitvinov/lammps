#ifndef LMP_UTILS_H
#define LMP_UTILS_H
namespace LAMMPS_NS {
class LAMMPS;
namespace utils {
void sfgets(const char *srcname, int srcline, char *s, int size, FILE *fp,
            const char *filename);
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
            bigint nmax, TYPE &nlo, TYPE &nhi);
char *strdup(const std::string &text);
std::string strip_style_suffix(const std::string &style, LAMMPS *lmp);
enum { NOCONVERT = 0, METAL2REAL = 1, REAL2METAL = 1 << 1 };
enum { UNKNOWN = 0, ENERGY };
} // namespace utils
} // namespace LAMMPS_NS
#endif
