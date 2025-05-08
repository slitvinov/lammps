#ifndef LMP_REGION_BLOCK_H
#define LMP_REGION_BLOCK_H
namespace LAMMPS_NS {
class RegBlock : public Region {
  friend class FixPour;

public:
  RegBlock(class LAMMPS *, int, char **);

protected:
  double xlo, xhi, ylo, yhi, zlo, zhi;
};
} // namespace LAMMPS_NS
#endif
