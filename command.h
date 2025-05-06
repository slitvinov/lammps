#ifndef LMP_COMMAND_H
#define LMP_COMMAND_H
namespace LAMMPS_NS {
class Command : protected Pointers {
public:
  Command(class LAMMPS *lmp) : Pointers(lmp){};
  virtual void command(int, char **) = 0;
};
} // namespace LAMMPS_NS
#endif
