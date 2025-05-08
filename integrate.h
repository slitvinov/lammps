#ifndef LMP_INTEGRATE_H
#define LMP_INTEGRATE_H
namespace LAMMPS_NS {
class Integrate : protected Pointers {
public:
  Integrate(class LAMMPS *, int, char **);
  virtual void init();
  virtual void setup(int flag) = 0;
  virtual void run(int) = 0;
  virtual void force_clear() = 0;
  virtual void cleanup() {}

protected:
  void ev_set(bigint);
};
} // namespace LAMMPS_NS
#endif
