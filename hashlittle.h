#ifndef LMP_HASHLITTLE_H
#define LMP_HASHLITTLE_H 
#include <cstddef>
#include <cstdint>
namespace LAMMPS_NS {
uint32_t hashlittle(const void *key, size_t length, uint32_t);
}
#endif
