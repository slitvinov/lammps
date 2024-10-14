#ifdef COMMAND_CLASS
CommandStyle(run, Run);
#else
#ifndef LMP_RUN_H
#define LMP_RUN_H
#include "command.h"
namespace LAMMPS_NS {
class Run : public Command {
public:
  Run(class LAMMPS *);
  void command(int, char **) override;
};
} // namespace LAMMPS_NS
#endif
#endif
