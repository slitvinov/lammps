#include "utils.h"
#include "arg_info.h"
#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "error.h"
#include "fix.h"
#include "fmt/chrono.h"
#include "input.h"
#include "label_map.h"
#include "memory.h"
#include "modify.h"
#include "text_file_reader.h"
#include "universe.h"
#include "update.h"
#include "variable.h"
#include <cctype>
#include <cerrno>
#include <cstring>
#include <ctime>
extern "C" {
static int re_match(const char *text, const char *pattern);
static int re_find(const char *text, const char *pattern, int *matchlen);
}
static void do_merge(int *idx, int *buf, int llo, int lhi, int rlo, int rhi, void *ptr,
                     int (*comp)(int, int, void *));
static void insertion_sort(int *index, int num, void *ptr, int (*comp)(int, int, void *));
using namespace LAMMPS_NS;
bool utils::strmatch(const std::string &text, const std::string &pattern)
{
  const int pos = re_match(text.c_str(), pattern.c_str());
  return (pos >= 0);
}
std::string utils::strfind(const std::string &text, const std::string &pattern)
{
  int matchlen;
  const int pos = re_find(text.c_str(), pattern.c_str(), &matchlen);
  if ((pos >= 0) && (matchlen > 0))
    return text.substr(pos, matchlen);
  else
    return "";
}
void utils::missing_cmd_args(const std::string &file, int line, const std::string &cmd,
                             Error *error)
{
  if (error) error->all(file, line, "Illegal {} command: missing argument(s)", cmd);
}
void utils::logmesg(LAMMPS *lmp, const std::string &mesg)
{
  if (lmp->screen) fputs(mesg.c_str(), lmp->screen);
  if (lmp->logfile) fputs(mesg.c_str(), lmp->logfile);
}
void utils::fmtargs_logmesg(LAMMPS *lmp, fmt::string_view format, fmt::format_args args)
{
  try {
    logmesg(lmp, fmt::vformat(format, args));
  } catch (fmt::format_error &e) {
    logmesg(lmp, std::string(e.what()) + "\n");
  }
}
std::string utils::errorurl(int errorcode)
{
  return fmt::format("\nFor more information see https://docs.lammps.org/err{:04d}", errorcode);
}
void utils::flush_buffers(LAMMPS *lmp)
{
  if (lmp->screen) fflush(lmp->screen);
  if (lmp->logfile) fflush(lmp->logfile);
  if (lmp->universe->uscreen) fflush(lmp->universe->uscreen);
  if (lmp->universe->ulogfile) fflush(lmp->universe->ulogfile);
}
std::string utils::getsyserror()
{
  return {strerror(errno)};
}
char *utils::fgets_trunc(char *buf, int size, FILE *fp)
{
  constexpr int MAXDUMMY = 256;
  char dummy[MAXDUMMY];
  char *ptr = fgets(buf, size, fp);
  if (!ptr) return nullptr;
  int n = strlen(buf);
  if (n < size - 1) {
    if (buf[n - 1] != '\n') {
      buf[n] = '\n';
      buf[n + 1] = '\0';
    }
    return buf;
  } else if (buf[n - 1] == '\n') {
    return buf;
  } else
    buf[size - 2] = '\n';
  do {
    ptr = fgets(dummy, MAXDUMMY, fp);
    if (ptr)
      n = strlen(ptr);
    else
      n = 0;
  } while (n == MAXDUMMY - 1 && ptr[MAXDUMMY - 1] != '\n');
  return buf;
}
void utils::sfgets(const char *srcname, int srcline, char *s, int size, FILE *fp,
                   const char *filename, Error *error)
{
  constexpr int MAXPATHLENBUF = 1024;
  char *rv = fgets(s, size, fp);
  if (rv == nullptr) {
    char buf[MAXPATHLENBUF];
    std::string errmsg;
    if (!filename) filename = platform::guesspath(fp, buf, MAXPATHLENBUF);
    if (feof(fp)) {
      errmsg = "Unexpected end of file while reading file '";
    } else if (ferror(fp)) {
      errmsg = "Unexpected error while reading file '";
    } else {
      errmsg = "Unexpected short read while reading file '";
    }
    errmsg += filename;
    errmsg += "'";
    if (error) error->one(srcname, srcline, errmsg);
    if (s) *s = '\0';
  }
}
void utils::sfread(const char *srcname, int srcline, void *s, size_t size, size_t num, FILE *fp,
                   const char *filename, Error *error)
{
  constexpr int MAXPATHLENBUF = 1024;
  size_t rv = fread(s, size, num, fp);
  if (rv != num) {
    char buf[MAXPATHLENBUF];
    std::string errmsg;
    if (!filename) filename = platform::guesspath(fp, buf, MAXPATHLENBUF);
    if (feof(fp)) {
      errmsg = "Unexpected end of file while reading file '";
    } else if (ferror(fp)) {
      errmsg = "Unexpected error while reading file '";
    } else {
      errmsg = "Unexpected short read while reading file '";
    }
    errmsg += filename;
    errmsg += "'";
    if (error) error->one(srcname, srcline, errmsg);
  }
}
int utils::read_lines_from_file(FILE *fp, int nlines, int nmax, char *buffer, int me, MPI_Comm comm)
{
  char *ptr = buffer;
  *ptr = '\0';
  if (me == 0) {
    if (fp) {
      for (int i = 0; i < nlines; i++) {
        ptr = fgets_trunc(ptr, nmax, fp);
        if (!ptr) break;
        ptr += strlen(ptr);
        *ptr = '\0';
      }
    }
  }
  int n = strlen(buffer);
  MPI_Bcast(&n, 1, MPI_INT, 0, comm);
  if (n == 0) return 1;
  MPI_Bcast(buffer, n + 1, MPI_CHAR, 0, comm);
  return 0;
}
int utils::logical(const char *file, int line, const std::string &str, bool do_abort, LAMMPS *lmp)
{
  if (str.empty()) {
    const char msg[] = "Expected boolean parameter instead of NULL or empty string "
                       "in input script or data file";
    if (do_abort)
      lmp->error->one(file, line, msg);
    else
      lmp->error->all(file, line, msg);
  }
  std::string buf(str);
  if (has_utf8(buf)) buf = utf8_subst(buf);
  int rv = 0;
  if ((buf == "yes") || (buf == "on") || (buf == "true") || (buf == "1")) {
    rv = 1;
  } else if ((buf == "no") || (buf == "off") || (buf == "false") || (buf == "0")) {
    rv = 0;
  } else {
    std::string msg("Expected boolean parameter instead of '");
    msg += buf + "' in input script or data file";
    if (do_abort)
      lmp->error->one(file, line, msg);
    else
      lmp->error->all(file, line, msg);
  }
  return rv;
}
int utils::logical(const char *file, int line, const char *str, bool do_abort, LAMMPS *lmp)
{
  if (str)
    return logical(file, line, std::string(str), do_abort, lmp);
  else
    return logical(file, line, std::string(""), do_abort, lmp);
}
double utils::numeric(const char *file, int line, const std::string &str, bool do_abort,
                      LAMMPS *lmp)
{
  if (str.empty()) {
    const char msg[] = "Expected floating point parameter instead of"
                       " NULL or empty string in input script or data file";
    if (do_abort)
      lmp->error->one(file, line, msg);
    else
      lmp->error->all(file, line, msg);
  }
  std::string buf(str);
  if (has_utf8(buf)) buf = utf8_subst(buf);
  if (!is_double(buf)) {
    std::string msg("Expected floating point parameter instead of '");
    msg += buf + "' in input script or data file";
    if (do_abort)
      lmp->error->one(file, line, msg);
    else
      lmp->error->all(file, line, msg);
  }
  return atof(buf.c_str());
}
double utils::numeric(const char *file, int line, const char *str, bool do_abort, LAMMPS *lmp)
{
  if (str)
    return numeric(file, line, std::string(str), do_abort, lmp);
  else
    return numeric(file, line, std::string(""), do_abort, lmp);
}
int utils::inumeric(const char *file, int line, const std::string &str, bool do_abort, LAMMPS *lmp)
{
  if (str.empty()) {
    const char msg[] = "Expected integer parameter instead of"
                       " NULL or empty string in input script or data file";
    if (do_abort)
      lmp->error->one(file, line, msg);
    else
      lmp->error->all(file, line, msg);
  }
  std::string buf(str);
  if (has_utf8(buf)) buf = utf8_subst(buf);
  if (!is_integer(buf)) {
    std::string msg("Expected integer parameter instead of '");
    msg += buf + "' in input script or data file";
    if (do_abort)
      lmp->error->one(file, line, msg);
    else
      lmp->error->all(file, line, msg);
  }
  return atoi(buf.c_str());
}
int utils::inumeric(const char *file, int line, const char *str, bool do_abort, LAMMPS *lmp)
{
  if (str)
    return inumeric(file, line, std::string(str), do_abort, lmp);
  else
    return inumeric(file, line, std::string(""), do_abort, lmp);
}
bigint utils::bnumeric(const char *file, int line, const std::string &str, bool do_abort,
                       LAMMPS *lmp)
{
  if (str.empty()) {
    const char msg[] = "Expected integer parameter instead of"
                       " NULL or empty string in input script or data file";
    if (do_abort)
      lmp->error->one(file, line, msg);
    else
      lmp->error->all(file, line, msg);
  }
  std::string buf(str);
  if (has_utf8(buf)) buf = utf8_subst(buf);
  if (!is_integer(buf)) {
    std::string msg("Expected integer parameter instead of '");
    msg += buf + "' in input script or data file";
    if (do_abort)
      lmp->error->one(file, line, msg);
    else
      lmp->error->all(file, line, msg);
  }
  return ATOBIGINT(buf.c_str());
}
bigint utils::bnumeric(const char *file, int line, const char *str, bool do_abort, LAMMPS *lmp)
{
  if (str)
    return bnumeric(file, line, std::string(str), do_abort, lmp);
  else
    return bnumeric(file, line, std::string(""), do_abort, lmp);
}
tagint utils::tnumeric(const char *file, int line, const std::string &str, bool do_abort,
                       LAMMPS *lmp)
{
  if (str.empty()) {
    const char msg[] = "Expected integer parameter instead of"
                       " NULL or empty string in input script or data file";
    if (do_abort)
      lmp->error->one(file, line, msg);
    else
      lmp->error->all(file, line, msg);
  }
  std::string buf(str);
  if (has_utf8(buf)) buf = utf8_subst(buf);
  if (!is_integer(buf)) {
    std::string msg("Expected integer parameter instead of '");
    msg += buf + "' in input script or data file";
    if (do_abort)
      lmp->error->one(file, line, msg);
    else
      lmp->error->all(file, line, msg);
  }
  return ATOTAGINT(buf.c_str());
}
tagint utils::tnumeric(const char *file, int line, const char *str, bool do_abort, LAMMPS *lmp)
{
  if (str)
    return tnumeric(file, line, std::string(str), do_abort, lmp);
  else
    return tnumeric(file, line, std::string(""), do_abort, lmp);
}
template <typename TYPE>
void utils::bounds(const char *file, int line, const std::string &str,
                   bigint nmin, bigint nmax, TYPE &nlo, TYPE &nhi, Error *error)
{
  nlo = nhi = -1;
  size_t found = str.find_first_not_of("*-0123456789");
  if (found != std::string::npos) {
    if (error) error->all(file, line, fmt::format("Invalid range string: {}", str));
    return;
  }
  found = str.find_first_of('*');
  if (found == std::string::npos) {
    nlo = nhi = strtol(str.c_str(), nullptr, 10);
  } else if (str.size() == 1) {
    nlo = nmin;
    nhi = nmax;
  } else if (found == 0) {
    nlo = nmin;
    nhi = strtol(str.substr(1).c_str(), nullptr, 10);
  } else if (str.size() == found + 1) {
    nlo = strtol(str.c_str(), nullptr, 10);
    nhi = nmax;
  } else {
    nlo = strtol(str.c_str(), nullptr, 10);
    nhi = strtol(str.substr(found + 1).c_str(), nullptr, 10);
  }
  if (error) {
    if ((nlo <= 0) || (nhi <= 0))
      error->all(file, line, fmt::format("Invalid range string: {}", str));
    if (nlo < nmin)
      error->all(file, line, fmt::format("Numeric index {} is out of bounds ({}-{})",
                                         nlo, nmin, nmax));
    else if (nhi > nmax)
      error->all(file, line, fmt::format("Numeric index {} is out of bounds ({}-{})",
                                         nhi, nmin, nmax));
    else if (nlo > nhi)
      error->all(file, line, fmt::format("Numeric index {} is out of bounds ({}-{})",
                                         nlo, nmin, nhi));
  }
}
template void utils::bounds<>(const char *, int, const std::string &,
                              bigint, bigint, int &, int &, Error *);
template void utils::bounds<>(const char *, int, const std::string &,
                              bigint, bigint, long &, long &, Error *);
template void utils::bounds<>(const char *, int, const std::string &,
                              bigint, bigint, long long &, long long &, Error *);
int utils::expand_args(const char *file, int line, int narg, char **arg, int mode, char **&earg,
                       LAMMPS *lmp)
{
  int iarg;
  char *ptr = nullptr;
  for (iarg = 0; iarg < narg; iarg++) {
    ptr = strchr(arg[iarg], '*');
    if (ptr) break;
  }
  if (!ptr) {
    earg = arg;
    return narg;
  }
  int maxarg = narg - iarg;
  earg = (char **) lmp->memory->smalloc(maxarg * sizeof(char *), "input:earg");
  int newarg = 0, expandflag, nlo, nhi, nmax;
  std::string id, wc, tail;
  for (iarg = 0; iarg < narg; iarg++) {
    std::string word(arg[iarg]);
    expandflag = 0;
    if (strmatch(word, "^[cf]_\\w+:\\w+:\\w+\\[\\d*\\*\\d*\\]")) {
      auto gridid = utils::parse_grid_id(FLERR, word, lmp->error);
      size_t first = gridid[2].find('[');
      size_t second = gridid[2].find(']', first + 1);
      id = gridid[2].substr(0, first);
      wc = gridid[2].substr(first + 1, second - first - 1);
      tail = gridid[2].substr(second + 1);
      if (gridid[0][0] == 'c') {
        auto compute = lmp->modify->get_compute_by_id(gridid[0].substr(2));
        if (compute && compute->pergrid_flag) {
          int dim = 0;
          int igrid = compute->get_grid_by_name(gridid[1], dim);
          if (igrid >= 0) {
            int ncol = 0;
            compute->get_griddata_by_name(igrid, id, ncol);
            nmax = ncol;
            expandflag = 1;
          }
        }
      } else if (gridid[0][0] == 'f') {
        auto fix = lmp->modify->get_fix_by_id(gridid[0].substr(2));
        if (fix && fix->pergrid_flag) {
          int dim = 0;
          int igrid = fix->get_grid_by_name(gridid[1], dim);
          if (igrid >= 0) {
            int ncol = 0;
            fix->get_griddata_by_name(igrid, id, ncol);
            nmax = ncol;
            expandflag = 1;
          }
        }
      }
      if (expandflag) {
        utils::bounds(file, line, wc, 1, nmax, nlo, nhi, lmp->error);
        if (newarg + nhi - nlo + 1 > maxarg) {
          maxarg += nhi - nlo + 1;
          earg = (char **) lmp->memory->srealloc(earg, maxarg * sizeof(char *), "input:earg");
        }
        for (int index = nlo; index <= nhi; index++) {
          earg[newarg] =
              utils::strdup(fmt::format("{}:{}:{}[{}]{}", gridid[0], gridid[1], id, index, tail));
          newarg++;
        }
      }
    } else if (strmatch(word, "^[cfv]_\\w+\\[\\d*\\*\\d*\\]") ||
               strmatch(word, "^[id]2_\\w+\\[\\d*\\*\\d*\\]")) {
      size_t first = word.find('[');
      size_t second = word.find(']', first + 1);
      if (word[1] == '2')
        id = word.substr(3, first - 3);
      else
        id = word.substr(2, first - 2);
      wc = word.substr(first + 1, second - first - 1);
      tail = word.substr(second + 1);
      if (word[0] == 'c') {
        auto compute = lmp->modify->get_compute_by_id(id);
        if (compute) {
          if (mode == 0 && compute->vector_flag) {
            nmax = compute->size_vector;
            expandflag = 1;
          } else if (mode == 1 && compute->array_flag) {
            nmax = compute->size_array_cols;
            expandflag = 1;
          } else if (compute->peratom_flag && compute->size_peratom_cols) {
            nmax = compute->size_peratom_cols;
            expandflag = 1;
          } else if (compute->local_flag && compute->size_local_cols) {
            nmax = compute->size_local_cols;
            expandflag = 1;
          }
        }
      } else if (word[0] == 'f') {
        auto fix = lmp->modify->get_fix_by_id(id);
        if (fix) {
          if (mode == 0 && fix->vector_flag) {
            nmax = fix->size_vector;
            expandflag = 1;
          } else if (mode == 1 && fix->array_flag) {
            nmax = fix->size_array_cols;
            expandflag = 1;
          } else if (fix->peratom_flag && fix->size_peratom_cols) {
            nmax = fix->size_peratom_cols;
            expandflag = 1;
          } else if (fix->local_flag && fix->size_local_cols) {
            nmax = fix->size_local_cols;
            expandflag = 1;
          }
        }
      } else if (word[0] == 'v') {
        int index = lmp->input->variable->find(id.c_str());
        if (index >= 0) {
          if (mode == 0 && lmp->input->variable->vectorstyle(index)) {
            utils::bounds(file, line, wc, 1, MAXSMALLINT, nlo, nhi, lmp->error);
            if (nhi < MAXSMALLINT) {
              nmax = nhi;
              expandflag = 1;
            }
          }
        }
      } else if ((word[0] == 'i') || (word[0] == 'd')) {
        int flag, cols;
        int icustom = lmp->atom->find_custom(id.c_str(), flag, cols);
        if ((icustom >= 0) && (mode == 1) && (cols > 0)) {
          if (((word[0] == 'i') && (flag == 0)) || ((word[0] == 'd') && (flag == 1))) {
            nmax = cols;
            expandflag = 1;
          }
        }
      }
      if (expandflag) {
        utils::bounds(file, line, wc, 1, nmax, nlo, nhi, lmp->error);
        if (newarg + nhi - nlo + 1 > maxarg) {
          maxarg += nhi - nlo + 1;
          earg = (char **) lmp->memory->srealloc(earg, maxarg * sizeof(char *), "input:earg");
        }
        for (int index = nlo; index <= nhi; index++) {
          if (word[1] == '2')
            earg[newarg] = utils::strdup(fmt::format("{}2_{}[{}]{}", word[0], id, index, tail));
          else
            earg[newarg] = utils::strdup(fmt::format("{}_{}[{}]{}", word[0], id, index, tail));
          newarg++;
        }
      }
    }
    if (!expandflag) {
      if (newarg == maxarg) {
        maxarg++;
        earg = (char **) lmp->memory->srealloc(earg, maxarg * sizeof(char *), "input:earg");
      }
      earg[newarg] = utils::strdup(word);
      newarg++;
    }
  }
  return newarg;
}
static const char *labeltypes[] = {"Atom", "Bond", "Angle", "Dihedral", "Improper"};
int utils::check_grid_reference(char *errstr, char *ref, int nevery,
                                char *&id, int &igrid, int &idata, int &index, LAMMPS *lmp)
{
  ArgInfo argi(ref, ArgInfo::COMPUTE | ArgInfo::FIX);
  index = argi.get_index1();
  auto name = argi.get_name();
  switch (argi.get_type()) {
    case ArgInfo::UNKNOWN: {
      lmp->error->all(FLERR,"%s grid reference %s is invalid",errstr,ref);
    } break;
    case ArgInfo::COMPUTE: {
      auto words = parse_grid_id(FLERR,name,lmp->error);
      const auto &idcompute = words[0];
      const auto &gname = words[1];
      const auto &dname = words[2];
      auto icompute = lmp->modify->get_compute_by_id(idcompute);
      if (!icompute) lmp->error->all(FLERR,"{} compute ID {} not found",errstr,idcompute);
      if (icompute->pergrid_flag == 0)
        lmp->error->all(FLERR,"{} compute {} does not compute per-grid info",errstr,idcompute);
      int dim;
      igrid = icompute->get_grid_by_name(gname,dim);
      if (igrid < 0)
        lmp->error->all(FLERR,"{} compute {} does not recognize grid name {}",errstr,idcompute,gname);
      int ncol;
      idata = icompute->get_griddata_by_name(igrid,dname,ncol);
      if (idata < 0)
        lmp->error->all(FLERR,"{} compute {} does not recognize data name {}",errstr,idcompute,dname);
      if (argi.get_dim() == 0 && ncol)
        lmp->error->all(FLERR,"{} compute {} data {} is not per-grid vector",errstr,idcompute,dname);
      if (argi.get_dim() && ncol == 0)
        lmp->error->all(FLERR,"{} compute {} data {} is not per-grid array",errstr,idcompute,dname);
      if (argi.get_dim() && argi.get_index1() > ncol)
        lmp->error->all(FLERR,"{} compute {} array {} is accessed out-of-range",errstr,idcompute,dname);
      id = utils::strdup(idcompute);
      return ArgInfo::COMPUTE;
    } break;
    case ArgInfo::FIX: {
      auto words = parse_grid_id(FLERR,name,lmp->error);
      const auto &idfix = words[0];
      const auto &gname = words[1];
      const auto &dname = words[2];
      auto ifix = lmp->modify->get_fix_by_id(idfix);
      if (!ifix) lmp->error->all(FLERR,"{} fix ID {} not found",errstr,idfix);
      if (ifix->pergrid_flag == 0)
        lmp->error->all(FLERR,"{} fix {} does not compute per-grid info",errstr,idfix);
      if (nevery % ifix->pergrid_freq)
        lmp->error->all(FLERR,"{} fix {} not computed at compatible time",errstr,idfix);
      int dim;
      igrid = ifix->get_grid_by_name(gname,dim);
      if (igrid < 0)
        lmp->error->all(FLERR,"{} fix {} does not recognize grid name {}",errstr,idfix,gname);
      int ncol;
      idata = ifix->get_griddata_by_name(igrid,dname,ncol);
      if (idata < 0)
        lmp->error->all(FLERR,"{} fix {} does not recognize data name {}",errstr,idfix,dname);
      if (argi.get_dim() == 0 && ncol)
        lmp->error->all(FLERR,"{} fix {} data {} is not per-grid vector",errstr,idfix,dname);
      if (argi.get_dim() > 0 && ncol == 0)
        lmp->error->all(FLERR,"{} fix {} data {} is not per-grid array",errstr,idfix,dname);
      if (argi.get_dim() > 0 && argi.get_index1() > ncol)
        lmp->error->all(FLERR,"{} fix {} array {} is accessed out-of-range",errstr,idfix,dname);
      id = utils::strdup(idfix);
      return ArgInfo::FIX;
    } break;
  }
  return -1;
}
std::vector<std::string> utils::parse_grid_id(const char *file, int line, const std::string &name,
                                              Error *error)
{
  auto words = Tokenizer(name, ":").as_vector();
  if (words.size() != 3) {
    if (error)
      error->all(file, line, "Grid ID {} must be 3 strings separated by 2 ':'characters", name);
    else
      return {"", "", ""};
  }
  return words;
}
char *utils::strdup(const std::string &text)
{
  auto tmp = new char[text.size() + 1];
  strcpy(tmp, text.c_str());
  return tmp;
}
std::string utils::lowercase(const std::string &text)
{
  std::string converted(text);
  for (auto &c : converted) c = ::tolower(c);
  return converted;
}
std::string utils::uppercase(const std::string &text)
{
  std::string converted(text);
  for (auto &c : converted) c = ::toupper(c);
  return converted;
}
std::string utils::trim(const std::string &line)
{
  int beg = re_match(line.c_str(), "\\S+");
  int end = re_match(line.c_str(), "\\s+$");
  if (beg < 0) beg = 0;
  if (end < 0) end = line.size();
  return line.substr(beg, end - beg);
}
std::string utils::trim_comment(const std::string &line)
{
  auto end = line.find('#');
  if (end != std::string::npos) { return line.substr(0, end); }
  return {line};
}
std::string utils::star_subst(const std::string &name, bigint step, int pad)
{
  auto star = name.find('*');
  if (star == std::string::npos) return name;
  return fmt::format("{}{:0{}}{}", name.substr(0, star), step, pad, name.substr(star + 1));
}
std::string utils::strip_style_suffix(const std::string &style, LAMMPS *lmp)
{
  std::string newstyle = style;
  return newstyle;
}
std::string utils::utf8_subst(const std::string &line)
{
  const auto *const in = (const unsigned char *) line.c_str();
  const int len = line.size();
  std::string out;
  for (int i = 0; i < len; ++i) {
    if ((in[i] & 0xe0U) == 0xc0U) {
      if ((i + 1) < len) {
        if ((in[i] == 0xc2U) && (in[i + 1] == 0xa0U)) out += ' ', ++i;
        if ((in[i] == 0xcbU) && (in[i + 1] == 0x96U)) out += '+', ++i;
        if ((in[i] == 0xcbU) && (in[i + 1] == 0x97U)) out += '-', ++i;
      }
    } else if ((in[i] & 0xf0U) == 0xe0U) {
      if ((i + 2) < len) {
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x80U)) out += ' ', i += 2;
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x81U)) out += ' ', i += 2;
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x82U)) out += ' ', i += 2;
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x83U)) out += ' ', i += 2;
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x84U)) out += ' ', i += 2;
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x85U)) out += ' ', i += 2;
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x86U)) out += ' ', i += 2;
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x87U)) out += ' ', i += 2;
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x88U)) out += ' ', i += 2;
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x89U)) out += ' ', i += 2;
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x8aU)) out += ' ', i += 2;
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x8bU)) out += ' ', i += 2;
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x98U)) out += '\'', i += 2;
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x99U)) out += '\'', i += 2;
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x9cU)) out += '"', i += 2;
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x9dU)) out += '"', i += 2;
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0xafU)) out += ' ', i += 2;
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x81U) && (in[i + 2] == 0xa0U)) out += ' ', i += 2;
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x81U) && (in[i + 2] == 0xa3U)) out += ' ', i += 2;
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x81U) && (in[i + 2] == 0xa4U)) out += '+', i += 2;
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x88U) && (in[i + 2] == 0x92U)) out += '-', i += 2;
        if ((in[i] == 0xefU) && (in[i + 1] == 0xbbU) && (in[i + 2] == 0xbfU)) out += ' ', i += 2;
      }
    } else if ((in[i] & 0xf8U) == 0xf0U) {
      if ((i + 3) < len) { ; }
    } else
      out += in[i];
  }
  return out;
}
size_t utils::count_words(const char *text)
{
  size_t count = 0;
  const char *buf = text;
  char c = *buf;
  while (c) {
    if (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\f') {
      c = *++buf;
      continue;
    };
    ++count;
    c = *++buf;
    while (c) {
      if (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\f') { break; }
      c = *++buf;
    }
  }
  return count;
}
size_t utils::count_words(const std::string &text)
{
  return utils::count_words(text.c_str());
}
size_t utils::count_words(const std::string &text, const std::string &separators)
{
  size_t count = 0;
  size_t start = text.find_first_not_of(separators);
  while (start != std::string::npos) {
    size_t end = text.find_first_of(separators, start);
    ++count;
    if (end == std::string::npos) {
      return count;
    } else {
      start = text.find_first_not_of(separators, end + 1);
    }
  }
  return count;
}
size_t utils::trim_and_count_words(const std::string &text, const std::string &separators)
{
  return utils::count_words(trim_comment(text), separators);
}
std::string utils::join_words(const std::vector<std::string> &words, const std::string &sep)
{
  std::string result;
  if (words.size() > 0) result = words[0];
  for (std::size_t i = 1; i < words.size(); ++i) result += sep + words[i];
  return result;
}
std::vector<std::string> utils::split_words(const std::string &text)
{
  std::vector<std::string> list;
  const char *buf = text.c_str();
  std::size_t beg = 0;
  std::size_t len = 0;
  std::size_t add = 0;
  char c = *buf;
  while (c) {
    if (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\f') {
      c = *++buf;
      ++beg;
      continue;
    };
    len = 0;
  quoted:
    if (c == '\'') {
      ++beg;
      add = 1;
      c = *++buf;
      while (((c != '\'') && (c != '\0')) || ((c == '\\') && (buf[1] == '\''))) {
        if ((c == '\\') && (buf[1] == '\'')) {
          ++buf;
          ++len;
        }
        c = *++buf;
        ++len;
      }
      if (c != '\'') ++len;
      c = *++buf;
    } else if ((c == '"') && (buf[1] == '"') && (buf[2] == '"') && (buf[3] != '"')) {
      len = 3;
      add = 1;
      buf += 3;
      c = *buf;
    } else if (c == '"') {
      ++beg;
      add = 1;
      c = *++buf;
      while (((c != '"') && (c != '\0')) || ((c == '\\') && (buf[1] == '"'))) {
        if ((c == '\\') && (buf[1] == '"')) {
          ++buf;
          ++len;
        }
        c = *++buf;
        ++len;
      }
      if (c != '"') ++len;
      c = *++buf;
    }
    while (true) {
      if ((c == '\'') || (c == '"')) goto quoted;
      if ((c == '\\') && ((buf[1] == '\'') || (buf[1] == '"'))) {
        ++buf;
        ++len;
        c = *++buf;
        ++len;
      }
      if ((c == ' ') || (c == '\t') || (c == '\r') || (c == '\n') || (c == '\f') || (c == '\0')) {
        list.push_back(text.substr(beg, len));
        beg += len + add;
        break;
      }
      c = *++buf;
      ++len;
    }
  }
  return list;
}
std::vector<std::string> utils::split_lines(const std::string &text)
{
  return Tokenizer(text, "\r\n").as_vector();
}
bool utils::is_integer(const std::string &str)
{
  if (str.empty()) return false;
  return strmatch(str, "^[+-]?\\d+$");
}
bool utils::is_double(const std::string &str)
{
  if (str.empty()) return false;
  return strmatch(str, "^[+-]?\\d+\\.?\\d*$") || strmatch(str, "^[+-]?\\d+\\.?\\d*[eE][+-]?\\d+$") ||
      strmatch(str, "^[+-]?\\d*\\.?\\d+$") || strmatch(str, "^[+-]?\\d*\\.?\\d+[eE][+-]?\\d+$");
}
bool utils::is_id(const std::string &str)
{
  if (str.empty()) return false;
  for (const auto &c : str) {
    if (isalnum(c) || (c == '_')) continue;
    return false;
  }
  return true;
}
int utils::is_type(const std::string &str)
{
  if (str.empty()) return -1;
  bool numeric = true;
  int nstar = 0;
  for (const auto &c : str) {
    if (isdigit(c)) continue;
    if (c == '*') {
      ++nstar;
      continue;
    }
    numeric = false;
  }
  if (numeric && (nstar < 2)) return 0;
  if (isdigit(str[0]) || (str[0] == '*') || (str[0] == '#')) return -1;
  if (str.find_first_of(" \t\r\n\f") != std::string::npos) return -1;
  if (has_utf8(utf8_subst(str))) return -1;
  return 1;
}
std::string utils::get_potential_file_path(const std::string &path)
{
  if (platform::file_is_readable(path)) {
    return path;
  } else {
    for (const auto &dir : platform::list_pathenv("LAMMPS_POTENTIALS")) {
      auto pot = platform::path_basename(path);
      auto filepath = platform::path_join(dir, pot);
      if (platform::file_is_readable(filepath)) return filepath;
    }
  }
  return "";
}
std::string utils::get_potential_date(const std::string &path, const std::string &potential_name)
{
  TextFileReader reader(path, potential_name);
  reader.ignore_comments = false;
  char *line = reader.next_line();
  if (line == nullptr) return "";
  Tokenizer words(line);
  while (words.has_next()) {
    if (words.next() == "DATE:") {
      if (words.has_next()) return words.next();
    }
  }
  return "";
}
std::string utils::get_potential_units(const std::string &path, const std::string &potential_name)
{
  TextFileReader reader(path, potential_name);
  reader.ignore_comments = false;
  char *line = reader.next_line();
  if (line == nullptr) return "";
  Tokenizer words(line);
  while (words.has_next()) {
    if (words.next() == "UNITS:") {
      if (words.has_next()) return words.next();
    }
  }
  return "";
}
int utils::get_supported_conversions(const int property)
{
  if (property == ENERGY)
    return METAL2REAL | REAL2METAL;
  else
    return NOCONVERT;
}
double utils::get_conversion_factor(const int property, const int conversion)
{
  if (property == ENERGY) {
    if (conversion == NOCONVERT) {
      return 1.0;
    } else if (conversion == METAL2REAL) {
      return 23.060549;
    } else if (conversion == REAL2METAL) {
      return 1.0 / 23.060549;
    }
  }
  return 0.0;
}
FILE *utils::open_potential(const std::string &name, LAMMPS *lmp, int *auto_convert)
{
  auto error = lmp->error;
  auto me = lmp->comm->me;
  std::string filepath = get_potential_file_path(name);
  if (!filepath.empty()) {
    std::string unit_style = lmp->update->unit_style;
    std::string date = get_potential_date(filepath, "potential");
    std::string units = get_potential_units(filepath, "potential");
    if (!date.empty() && (me == 0))
      logmesg(lmp, "Reading potential file {} with DATE: {}\n", name, date);
    if (auto_convert == nullptr) {
      if (!units.empty() && (units != unit_style) && (me == 0)) {
        error->one(FLERR, "Potential file {} requires {} units but {} units are in use", name,
                   units, unit_style);
        return nullptr;
      }
    } else {
      if (units.empty() || units == unit_style) {
        *auto_convert = NOCONVERT;
      } else {
        if ((units == "metal") && (unit_style == "real") && (*auto_convert & METAL2REAL)) {
          *auto_convert = METAL2REAL;
        } else if ((units == "real") && (unit_style == "metal") && (*auto_convert & REAL2METAL)) {
          *auto_convert = REAL2METAL;
        } else {
          error->one(FLERR, "Potential file {} requires {} units but {} units are in use", name,
                     units, unit_style);
          return nullptr;
        }
      }
      if ((*auto_convert != NOCONVERT) && (me == 0))
        error->warning(FLERR, "Converting potential file in {} units to {} units", units,
                       unit_style);
    }
    return fopen(filepath.c_str(), "r");
  }
  return nullptr;
}
double utils::timespec2seconds(const std::string &timespec)
{
  double vals[3];
  int i = 0;
  if (timespec == "off") return -1.0;
  if (timespec == "unlimited") return -1.0;
  vals[0] = vals[1] = vals[2] = 0;
  ValueTokenizer values(timespec, ":");
  try {
    for (i = 0; i < 3; i++) {
      if (!values.has_next()) break;
      vals[i] = values.next_int();
    }
  } catch (TokenizerException &) {
    return -1.0;
  }
  if (i == 3)
    return (vals[0] * 60 + vals[1]) * 60 + vals[2];
  else if (i == 2)
    return vals[0] * 60 + vals[1];
  return vals[0];
}
int utils::date2num(const std::string &date)
{
  std::size_t found = date.find_first_not_of("0123456789 ");
  int num = strtol(date.substr(0, found).c_str(), nullptr, 10);
  auto month = date.substr(found);
  found = month.find_first_of("0123456789 ");
  num += strtol(month.substr(found).c_str(), nullptr, 10) * 10000;
  if (num < 1000000) num += 20000000;
  if (strmatch(month, "^Jan"))
    num += 100;
  else if (strmatch(month, "^Feb"))
    num += 200;
  else if (strmatch(month, "^Mar"))
    num += 300;
  else if (strmatch(month, "^Apr"))
    num += 400;
  else if (strmatch(month, "^May"))
    num += 500;
  else if (strmatch(month, "^Jun"))
    num += 600;
  else if (strmatch(month, "^Jul"))
    num += 700;
  else if (strmatch(month, "^Aug"))
    num += 800;
  else if (strmatch(month, "^Sep"))
    num += 900;
  else if (strmatch(month, "^Oct"))
    num += 1000;
  else if (strmatch(month, "^Nov"))
    num += 1100;
  else if (strmatch(month, "^Dec"))
    num += 1200;
  return num;
}
std::string utils::current_date()
{
  time_t tv = time(nullptr);
  std::tm today = fmt::localtime(tv);
  return fmt::format("{:%Y-%m-%d}", today);
}
int utils::binary_search(const double needle, const int n, const double *haystack)
{
  int lo = 0;
  int hi = n - 1;
  if (needle < haystack[lo]) return lo;
  if (needle >= haystack[hi]) return hi;
  int index = (lo + hi) / 2;
  while (lo < hi - 1) {
    if (needle < haystack[index])
      hi = index;
    else if (needle >= haystack[index])
      lo = index;
    index = (lo + hi) / 2;
  }
  return index;
}
void utils::merge_sort(int *index, int num, void *ptr, int (*comp)(int, int, void *))
{
  if (num < 2) return;
  int chunk, i, j;
  chunk = 64;
  for (i = 0; i < num; i += chunk) {
    j = (i + chunk > num) ? num - i : chunk;
    insertion_sort(index + i, j, ptr, comp);
  }
  if (chunk >= num) return;
  int *buf = new int[num];
  int *dest = index;
  int *hold = buf;
  while (chunk < num) {
    int m;
    int *tmp = dest;
    dest = hold;
    hold = tmp;
    for (i = 0; i < num - 1; i += 2 * chunk) {
      j = i + 2 * chunk;
      if (j > num) j = num;
      m = i + chunk;
      if (m > num) m = num;
      do_merge(dest, hold, i, m, m, j, ptr, comp);
    }
    for (; i < num; i++) dest[i] = hold[i];
    chunk *= 2;
  }
  if (dest == buf) memcpy(index, buf, sizeof(int) * num);
  delete[] buf;
}
void insertion_sort(int *index, int num, void *ptr, int (*comp)(int, int, void *))
{
  if (num < 2) return;
  for (int i = 1; i < num; ++i) {
    int tmp = index[i];
    for (int j = i - 1; j >= 0; --j) {
      if ((*comp)(index[j], tmp, ptr) > 0) {
        index[j + 1] = index[j];
      } else {
        index[j + 1] = tmp;
        break;
      }
      if (j == 0) index[0] = tmp;
    }
  }
}
static void do_merge(int *idx, int *buf, int llo, int lhi, int rlo, int rhi, void *ptr,
                     int (*comp)(int, int, void *))
{
  int i = llo;
  int l = llo;
  int r = rlo;
  while ((l < lhi) && (r < rhi)) {
    if ((*comp)(buf[l], buf[r], ptr) < 0)
      idx[i++] = buf[l++];
    else
      idx[i++] = buf[r++];
  }
  while (l < lhi) idx[i++] = buf[l++];
  while (r < rhi) idx[i++] = buf[r++];
}
extern "C" {
typedef struct regex_t *re_t;
typedef struct regex_context_t *re_ctx_t;
static re_t re_compile(re_ctx_t context, const char *pattern);
static int re_matchp(const char *text, re_t pattern, int *matchlen);
#define MAX_REGEXP_OBJECTS 256
#define MAX_CHAR_CLASS_LEN 256
enum {
  RX_UNUSED,
  RX_DOT,
  RX_BEGIN,
  RX_END,
  RX_QUESTIONMARK,
  RX_STAR,
  RX_PLUS,
  RX_CHAR,
  RX_CHAR_CLASS,
  RX_INV_CHAR_CLASS,
  RX_DIGIT,
  RX_NOT_DIGIT,
  RX_INTEGER,
  RX_NOT_INTEGER,
  RX_FLOAT,
  RX_NOT_FLOAT,
  RX_ALPHA,
  RX_NOT_ALPHA,
  RX_WHITESPACE,
  RX_NOT_WHITESPACE
};
typedef struct regex_t {
  unsigned char type;
  union {
    unsigned char ch;
    unsigned char *ccl;
  } u;
} regex_t;
typedef struct regex_context_t {
  regex_t re_compiled[MAX_REGEXP_OBJECTS];
  unsigned char ccl_buf[MAX_CHAR_CLASS_LEN];
} regex_context_t;
int re_match(const char *text, const char *pattern)
{
  regex_context_t context;
  int dummy;
  return re_matchp(text, re_compile(&context, pattern), &dummy);
}
int re_find(const char *text, const char *pattern, int *matchlen)
{
  regex_context_t context;
  return re_matchp(text, re_compile(&context, pattern), matchlen);
}
static int matchpattern(regex_t *pattern, const char *text, int *matchlen);
static int matchcharclass(char c, const char *str);
static int matchstar(regex_t p, regex_t *pattern, const char *text, int *matchlen);
static int matchplus(regex_t p, regex_t *pattern, const char *text, int *matchlen);
static int matchone(regex_t p, char c);
static int matchdigit(char c);
static int matchint(char c);
static int matchfloat(char c);
static int matchalpha(char c);
static int matchwhitespace(char c);
static int matchmetachar(char c, const char *str);
static int matchrange(char c, const char *str);
static int matchdot(char c);
static int ismetachar(char c);
int re_matchp(const char *text, re_t pattern, int *matchlen)
{
  *matchlen = 0;
  if (pattern != nullptr) {
    if (pattern[0].type == RX_BEGIN) {
      return ((matchpattern(&pattern[1], text, matchlen)) ? 0 : -1);
    } else {
      int idx = -1;
      do {
        idx += 1;
        if (matchpattern(pattern, text, matchlen)) {
          if (text[0] == '\0') return -1;
          return idx;
        }
      } while (*text++ != '\0');
    }
  }
  return -1;
}
re_t re_compile(re_ctx_t context, const char *pattern)
{
  regex_t *const re_compiled = context->re_compiled;
  unsigned char *const ccl_buf = context->ccl_buf;
  int ccl_bufidx = 1;
  char c;
  int i = 0;
  int j = 0;
  while (pattern[i] != '\0' && (j + 1 < MAX_REGEXP_OBJECTS)) {
    c = pattern[i];
    switch (c) {
      case '^': {
        re_compiled[j].type = RX_BEGIN;
      } break;
      case '$': {
        re_compiled[j].type = RX_END;
      } break;
      case '.': {
        re_compiled[j].type = RX_DOT;
      } break;
      case '*': {
        re_compiled[j].type = RX_STAR;
      } break;
      case '+': {
        re_compiled[j].type = RX_PLUS;
      } break;
      case '?': {
        re_compiled[j].type = RX_QUESTIONMARK;
      } break;
      case '\\': {
        if (pattern[i + 1] != '\0') {
          i += 1;
          switch (pattern[i]) {
            case 'd': {
              re_compiled[j].type = RX_DIGIT;
            } break;
            case 'D': {
              re_compiled[j].type = RX_NOT_DIGIT;
            } break;
            case 'i': {
              re_compiled[j].type = RX_INTEGER;
            } break;
            case 'I': {
              re_compiled[j].type = RX_NOT_INTEGER;
            } break;
            case 'f': {
              re_compiled[j].type = RX_FLOAT;
            } break;
            case 'F': {
              re_compiled[j].type = RX_NOT_FLOAT;
            } break;
            case 'w': {
              re_compiled[j].type = RX_ALPHA;
            } break;
            case 'W': {
              re_compiled[j].type = RX_NOT_ALPHA;
            } break;
            case 's': {
              re_compiled[j].type = RX_WHITESPACE;
            } break;
            case 'S': {
              re_compiled[j].type = RX_NOT_WHITESPACE;
            } break;
            default: {
              re_compiled[j].type = RX_CHAR;
              re_compiled[j].u.ch = pattern[i];
            } break;
          }
        }
      } break;
      case '[': {
        int buf_begin = ccl_bufidx;
        if (pattern[i + 1] == '^') {
          re_compiled[j].type = RX_INV_CHAR_CLASS;
          i += 1;
          if (pattern[i + 1] == 0)
          {
            return nullptr;
          }
        } else {
          re_compiled[j].type = RX_CHAR_CLASS;
        }
        while ((pattern[++i] != ']') && (pattern[i] != '\0')) {
          if (pattern[i] == '\\') {
            if (ccl_bufidx >= MAX_CHAR_CLASS_LEN - 1) { return nullptr; }
            if (pattern[i + 1] == 0)
            {
              return nullptr;
            }
            ccl_buf[ccl_bufidx++] = pattern[i++];
          } else if (ccl_bufidx >= MAX_CHAR_CLASS_LEN) {
            return nullptr;
          }
          ccl_buf[ccl_bufidx++] = pattern[i];
        }
        if (ccl_bufidx >= MAX_CHAR_CLASS_LEN) {
          return nullptr;
        }
        ccl_buf[ccl_bufidx++] = 0;
        re_compiled[j].u.ccl = &ccl_buf[buf_begin];
      } break;
      default: {
        re_compiled[j].type = RX_CHAR;
        re_compiled[j].u.ch = c;
      } break;
    }
    if (pattern[i] == 0) { return nullptr; }
    i += 1;
    j += 1;
  }
  re_compiled[j].type = RX_UNUSED;
  return (re_t) re_compiled;
}
static int matchdigit(char c)
{
  return isdigit(c);
}
static int matchint(char c)
{
  return (matchdigit(c) || (c == '-') || (c == '+'));
}
static int matchfloat(char c)
{
  return (matchint(c) || (c == '.') || (c == 'e') || (c == 'E'));
}
static int matchalpha(char c)
{
  return isalpha(c);
}
static int matchwhitespace(char c)
{
  return isspace(c);
}
static int matchalphanum(char c)
{
  return ((c == '_') || matchalpha(c) || matchdigit(c));
}
static int matchrange(char c, const char *str)
{
  return ((c != '-') && (str[0] != '\0') && (str[0] != '-') && (str[1] == '-') &&
          (str[1] != '\0') && (str[2] != '\0') && ((c >= str[0]) && (c <= str[2])));
}
static int matchdot(char c)
{
#if defined(RE_DOT_MATCHES_NEWLINE) && (RE_DOT_MATCHES_NEWLINE == 1)
  (void) c;
  return 1;
#else
  return c != '\n' && c != '\r';
#endif
}
static int ismetachar(char c)
{
  return ((c == 's') || (c == 'S') || (c == 'w') || (c == 'W') || (c == 'd') || (c == 'D'));
}
static int matchmetachar(char c, const char *str)
{
  switch (str[0]) {
    case 'd':
      return matchdigit(c);
    case 'D':
      return !matchdigit(c);
    case 'i':
      return matchint(c);
    case 'I':
      return !matchint(c);
    case 'f':
      return matchfloat(c);
    case 'F':
      return !matchfloat(c);
    case 'w':
      return matchalphanum(c);
    case 'W':
      return !matchalphanum(c);
    case 's':
      return matchwhitespace(c);
    case 'S':
      return !matchwhitespace(c);
    default:
      return (c == str[0]);
  }
}
static int matchcharclass(char c, const char *str)
{
  do {
    if (matchrange(c, str)) {
      return 1;
    } else if (str[0] == '\\') {
      str += 1;
      if (matchmetachar(c, str)) {
        return 1;
      } else if ((c == str[0]) && !ismetachar(c)) {
        return 1;
      }
    } else if (c == str[0]) {
      if (c == '-') {
        return ((str[-1] == '\0') || (str[1] == '\0'));
      } else {
        return 1;
      }
    }
  } while (*str++ != '\0');
  return 0;
}
static int matchone(regex_t p, char c)
{
  switch (p.type) {
    case RX_DOT:
      return matchdot(c);
    case RX_CHAR_CLASS:
      return matchcharclass(c, (const char *) p.u.ccl);
    case RX_INV_CHAR_CLASS:
      return !matchcharclass(c, (const char *) p.u.ccl);
    case RX_DIGIT:
      return matchdigit(c);
    case RX_NOT_DIGIT:
      return !matchdigit(c);
    case RX_INTEGER:
      return matchint(c);
    case RX_NOT_INTEGER:
      return !matchint(c);
    case RX_FLOAT:
      return matchfloat(c);
    case RX_NOT_FLOAT:
      return !matchfloat(c);
    case RX_ALPHA:
      return matchalphanum(c);
    case RX_NOT_ALPHA:
      return !matchalphanum(c);
    case RX_WHITESPACE:
      return matchwhitespace(c);
    case RX_NOT_WHITESPACE:
      return !matchwhitespace(c);
    default:
      return (p.u.ch == c);
  }
}
static int matchstar(regex_t p, regex_t *pattern, const char *text, int *matchlen)
{
  int prelen = *matchlen;
  const char *prepos = text;
  while ((text[0] != '\0') && matchone(p, *text)) {
    text++;
    (*matchlen)++;
  }
  while (text >= prepos) {
    if (matchpattern(pattern, text--, matchlen)) return 1;
    (*matchlen)--;
  }
  *matchlen = prelen;
  return 0;
}
static int matchplus(regex_t p, regex_t *pattern, const char *text, int *matchlen)
{
  const char *prepos = text;
  while ((text[0] != '\0') && matchone(p, *text)) {
    text++;
    (*matchlen)++;
  }
  while (text > prepos) {
    if (matchpattern(pattern, text--, matchlen)) return 1;
    (*matchlen)--;
  }
  return 0;
}
static int matchquestion(regex_t p, regex_t *pattern, const char *text, int *matchlen)
{
  if (p.type == RX_UNUSED) return 1;
  if (matchpattern(pattern, text, matchlen)) return 1;
  if (*text && matchone(p, *text++)) {
    if (matchpattern(pattern, text, matchlen)) {
      (*matchlen)++;
      return 1;
    }
  }
  return 0;
}
static int matchpattern(regex_t *pattern, const char *text, int *matchlen)
{
  int pre = *matchlen;
  do {
    if ((pattern[0].type == RX_UNUSED) || (pattern[1].type == RX_QUESTIONMARK)) {
      return matchquestion(pattern[0], &pattern[2], text, matchlen);
    } else if (pattern[1].type == RX_STAR) {
      return matchstar(pattern[0], &pattern[2], text, matchlen);
    } else if (pattern[1].type == RX_PLUS) {
      return matchplus(pattern[0], &pattern[2], text, matchlen);
    } else if ((pattern[0].type == RX_END) && pattern[1].type == RX_UNUSED) {
      return (text[0] == '\0');
    }
    (*matchlen)++;
  } while ((text[0] != '\0') && matchone(*pattern++, *text++));
  *matchlen = pre;
  return 0;
}
}
