#ifndef LMP_CREATE_BOX_H
#define LMP_CREATE_BOX_H
namespace LAMMPS_NS {
class CreateBox : public Command {
public:
  CreateBox(class LAMMPS *);
  void command(int, char **) override;
};
} // namespace LAMMPS_NS
#endif
