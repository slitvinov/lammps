#ifndef LMP_VERLET_H
#define LMP_VERLET_H
namespace LAMMPS_NS {
class Verlet : protected Pointers {
public:
  Verlet(class LAMMPS *, int, char **);
  void init();
  void setup(int flag);
  void run(int);
  void force_clear();
  void cleanup();

protected:
  int torqueflag, extraflag;
};
} // namespace LAMMPS_NS
#endif
