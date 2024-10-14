#ifdef COMMAND_CLASS
CommandStyle(create_box, CreateBox);
#else
#ifndef LMP_CREATE_BOX_H
#define LMP_CREATE_BOX_H
#include "command.h"
namespace LAMMPS_NS {
class CreateBox : public Command {
public:
  CreateBox(class LAMMPS *);
  void command(int, char **) override;
};
} // namespace LAMMPS_NS
#endif
#endif
