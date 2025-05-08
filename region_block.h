#ifndef LMP_REGION_BLOCK_H
#define LMP_REGION_BLOCK_H
namespace LAMMPS_NS {
class RegBlock {
public:
  RegBlock(class LAMMPS *, int, char **);
  double xlo, xhi, ylo, yhi, zlo, zhi;
};
} // namespace LAMMPS_NS
#endif
