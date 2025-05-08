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
};
} // namespace LAMMPS_NS
#endif
