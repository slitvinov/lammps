#ifdef COMMAND_CLASS
CommandStyle(set,Set);
#else
#ifndef LMP_SET_H
#define LMP_SET_H 
#include "command.h"
namespace LAMMPS_NS {
class Set : public Command {
 public:
  Set(class LAMMPS *lmp) : Command(lmp){};
  void command(int, char **) override;
 private:
  char *id;
  int *select;
  int style, ivalue, newtype, count, index_custom, icol_custom;
  int ximage, yimage, zimage, ximageflag, yimageflag, zimageflag;
  int cc_index;
  bigint nsubset;
  double dvalue, xvalue, yvalue, zvalue, wvalue, fraction;
  int varflag, varflag1, varflag2, varflag3, varflag4;
  int ivar1, ivar2, ivar3, ivar4;
  double *vec1, *vec2, *vec3, *vec4;
  int discflag;
  void selection(int);
  void set(int);
  void varparse(const char *, int);
};
}
#endif
#endif
