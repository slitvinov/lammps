#include "random_park.h"
#include <cmath>
using namespace LAMMPS_NS;
#define IA 16807
#define IM 2147483647
#define AM (1.0 / IM)
#define IQ 127773
#define IR 2836
RanPark::RanPark(LAMMPS *lmp, int seed_init) : Pointers(lmp) {
  seed = seed_init;
  save = 0;
}
double RanPark::uniform() {
  int k = seed / IQ;
  seed = IA * (seed - k * IQ) - IR * k;
  if (seed < 0)
    seed += IM;
  double ans = AM * seed;
  return ans;
}
