#ifdef COMMAND_CLASS
CommandStyle(create_atoms,CreateAtoms);
#else
#ifndef LMP_CREATE_ATOMS_H
#define LMP_CREATE_ATOMS_H 
#include "command.h"
namespace LAMMPS_NS {
class CreateAtoms : public Command {
 public:
  CreateAtoms(class LAMMPS *);
  void command(int, char **) override;
 private:
  int ntype, style, mode, nbasis, seed;
  bigint nrandom;
  int remapflag;
  int maxtry;
  int quat_user;
  int overlapflag;
  double overlap;
  int subsetflag;
  bigint nsubset;
  double subsetfrac;
  int *basistype;
  double xone[3], quatone[4];
  double radthresh, radscale, mesh_density;
  int varflag, vvar, xvar, yvar, zvar;
  char *vstr, *xstr, *ystr, *zstr;
  char *xstr_copy, *ystr_copy, *zstr_copy;
  int ilo, ihi, jlo, jhi, klo, khi;
  int nlatt;
  int nlatt_overflow;
  int *flag;
  int *next;
  int mesh_style;
  class Region *region;
  class RanMars *ranmol;
  class RanMars *ranlatt;
  int triclinic;
  double sublo[3], subhi[3];
  void add_single();
  void add_random();
  void loop_lattice(int);
  int vartest(double *);
};
}
#endif
#endif
