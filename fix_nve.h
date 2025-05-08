#ifndef LMP_FIX_NVE_H
#define LMP_FIX_NVE_H
namespace LAMMPS_NS {
class FixNVE : public Fix {
public:
  FixNVE(class LAMMPS *, int, char **);
  void initial_integrate(int) override;
  void final_integrate() override;
};
} // namespace LAMMPS_NS
#endif
