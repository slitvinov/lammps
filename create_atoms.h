#ifndef LMP_CREATE_ATOMS_H
#define LMP_CREATE_ATOMS_H
namespace LAMMPS_NS {
class CreateAtoms : public Command {
public:
  CreateAtoms(class LAMMPS *);
  void command(int, char **) override;

private:
  int ntype, seed;
  bigint nrandom;
  double overlap;
  bigint nsubset;
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
