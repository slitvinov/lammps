#include "arg_info.h"
#include <cstring>
#include <stdexcept>
using namespace LAMMPS_NS;
ArgInfo::ArgInfo(const std::string &arg, int allowed)
    : type(NONE), dim(0), index1(-1), index2(-1) {
  if (((arg.size() > 3) && (arg[1] == '2') && (arg[2] == '_')) ||
      ((arg.size() > 2) && (arg[1] == '_'))) {
    if ((arg[0] == 'c') && (allowed & COMPUTE))
      type = COMPUTE;
    else if ((arg[0] == 'f') && (allowed & FIX))
      type = FIX;
    else if ((arg[0] == 'v') && (allowed & VARIABLE))
      type = VARIABLE;
    else if ((arg[0] == 'd') && (allowed & DNAME))
      type = DNAME;
    else if ((arg[0] == 'i') && (allowed & INAME))
      type = INAME;
    else {
      index1 = 0;
      name = arg;
      return;
    }
    const int offset = (arg[1] == '_') ? 2 : 3;
    std::size_t has_idx1 = arg.find('[', offset);
    if (has_idx1 != std::string::npos) {
      name = arg.substr(offset, has_idx1 - offset);
      dim = 1;
      std::size_t has_idx2 = arg.find('[', has_idx1 + 1);
      if (has_idx2 != std::string::npos) {
        dim = 2;
        if (arg[arg.size() - 1] != ']') {
          type = UNKNOWN;
        } else {
          try {
            index2 = std::stoi(
                arg.substr(has_idx2 + 1, arg.size() - (has_idx2 + 2)));
          } catch (std::invalid_argument &) {
            type = UNKNOWN;
          }
        }
      } else
        has_idx2 = arg.size();
      if ((arg[has_idx2 - 1] != ']') ||
          ((dim == 1) && (arg.find(']') != has_idx2 - 1))) {
        type = UNKNOWN;
      } else {
        try {
          index1 =
              std::stoi(arg.substr(has_idx1 + 1, arg.size() - (has_idx1 + 2)));
        } catch (std::invalid_argument &) {
          type = UNKNOWN;
        }
      }
    } else {
      index1 = 0;
      name = arg.substr(offset);
    }
  } else {
    index1 = 0;
    name = arg;
  }
}
char *ArgInfo::copy_name() {
  auto dest = new char[name.size() + 1];
  strcpy(dest, name.c_str());
  return dest;
}
