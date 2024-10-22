#ifndef LMP_ERROR_H
#define LMP_ERROR_H
#include "pointers.h"
namespace LAMMPS_NS {
class Error : protected Pointers {
public:
  Error(class LAMMPS *);
  [[noreturn]] void universe_all(const std::string &, int, const std::string &);
  [[noreturn]] void universe_one(const std::string &, int, const std::string &);
  void universe_warn(const std::string &, int, const std::string &);
  [[noreturn]] void all(const std::string &, int, const std::string &);
  template <typename... Args>
  void all(const std::string &file, int line, const std::string &format,
           Args &&... args) {
    _all(file, line, format, fmt::make_format_args(args...));
  }
  [[noreturn]] void one(const std::string &, int, const std::string &);
  template <typename... Args>
  void one(const std::string &file, int line, const std::string &format,
           Args &&... args) {
    _one(file, line, format, fmt::make_format_args(args...));
  }
  void warning(const std::string &, int, const std::string &);
  template <typename... Args>
  void warning(const std::string &file, int line, const std::string &format,
               Args &&... args) {
    _warning(file, line, format, fmt::make_format_args(args...));
  }
  void message(const std::string &, int, const std::string &);
  template <typename... Args>
  void message(const std::string &file, int line, const std::string &format,
               Args &&... args) {
    _message(file, line, format, fmt::make_format_args(args...));
  }
  [[noreturn]] void done(int = 0);
  int get_numwarn() const { return numwarn; }
  int get_maxwarn() const { return maxwarn; }
  void set_numwarn(int val) { numwarn = val; }
  void set_maxwarn(int val) { maxwarn = val; }
  void set_allwarn(int val) { allwarn = val; }

private:
  int numwarn, maxwarn, allwarn;
  [[noreturn]] void _all(const std::string &, int, fmt::string_view,
                         fmt::format_args args);
  [[noreturn]] void _one(const std::string &, int, fmt::string_view,
                         fmt::format_args args);
  void _warning(const std::string &, int, fmt::string_view,
                fmt::format_args args);
  void _message(const std::string &, int, fmt::string_view,
                fmt::format_args args);
};
} // namespace LAMMPS_NS
#endif
