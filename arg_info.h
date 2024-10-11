#ifndef LMP_ARG_INFO_H
#define LMP_ARG_INFO_H 
#include <string>
namespace LAMMPS_NS {
class ArgInfo {
 public:
  enum ArgTypes {
    ERROR =-2,
    UNKNOWN =-1,
    NONE = 0,
    X = 1<<0,
    V = 1<<1,
    F = 1<<2,
    COMPUTE = 1<<3,
    FIX = 1<<4,
    VARIABLE = 1<<5,
    KEYWORD = 1<<6,
    TYPE = 1<<7,
    MOLECULE = 1<<8,
    DNAME = 1<<9,
    INAME = 1<<10,
    DENSITY_NUMBER = 1<<11,
    DENSITY_MASS = 1<<12,
    MASS = 1<<13,
    TEMPERATURE = 1<<14,
    BIN1D = 1<<15,
    BIN2D = 1<<16,
    BIN3D = 1<<17,
    BINSPHERE = 1<<18,
    BINCYLINDER = 1<<19
  };
  ArgInfo(const std::string &arg, int allowed = COMPUTE | FIX | VARIABLE);
  virtual ~ArgInfo() {}
 public:
  int get_type() const { return type; }
  int get_dim() const { return dim; }
  int get_index1() const { return index1; }
  int get_index2() const { return index2; }
  const char *get_name() const { return name.c_str(); }
  char *copy_name();
 private:
  std::string name;
  int type, dim, index1, index2;
  ArgInfo() {}
  ArgInfo(const ArgInfo &) {}
  void operator=(const ArgInfo &) {}
};
}
#endif
