#ifndef LMP_CREATE_ATOMS_H
#define LMP_CREATE_ATOMS_H
namespace LAMMPS_NS {
class CreateAtoms : public Command {
public:
  CreateAtoms(class LAMMPS *);
  void command(int, char **) override;

private:
  int ntype, style, seed;
  bigint nrandom;
  int quat_user;
  double overlap;
  bigint nsubset;
  double xone[3], quatone[4];
  int varflag, vvar, xvar, yvar, zvar;
  int ilo, ihi, jlo, jhi, klo, khi;
  int *flag;
  int *next;
  class Region *region;
  double sublo[3], subhi[3];
  void add_random();
};
} // namespace LAMMPS_NS
#endif
