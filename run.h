#ifndef LMP_RUN_H
#define LMP_RUN_H
namespace LAMMPS_NS {
class Run : public Command {
public:
  Run(class LAMMPS *);
  void command(int, char **) override;
};
} // namespace LAMMPS_NS
#endif
