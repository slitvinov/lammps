#include "random_park.h"
#include "error.h"
#include <cmath>
using namespace LAMMPS_NS;
#define IA 16807
#define IM 2147483647
#define AM (1.0 / IM)
#define IQ 127773
#define IR 2836
RanPark::RanPark(LAMMPS *lmp, int seed_init) : Pointers(lmp)
{
  if (seed_init <= 0) error->one(FLERR, "Invalid seed for Park random # generator");
  seed = seed_init;
  save = 0;
}
double RanPark::uniform()
{
  int k = seed / IQ;
  seed = IA * (seed - k * IQ) - IR * k;
  if (seed < 0) seed += IM;
  double ans = AM * seed;
  return ans;
}
double RanPark::gaussian()
{
  double first, v1, v2, rsq, fac;
  if (!save) {
    do {
      v1 = 2.0 * uniform() - 1.0;
      v2 = 2.0 * uniform() - 1.0;
      rsq = v1 * v1 + v2 * v2;
    } while ((rsq >= 1.0) || (rsq == 0.0));
    fac = sqrt(-2.0 * log(rsq) / rsq);
    second = v1 * fac;
    first = v2 * fac;
    save = 1;
  } else {
    first = second;
    save = 0;
  }
  return first;
}
void RanPark::reset(int seed_init)
{
  if (seed_init <= 0) error->all(FLERR, "Invalid seed for Park random # generator");
  seed = seed_init;
  save = 0;
}
void RanPark::reset(int ibase, double *coord)
{
  int i;
  auto str = (char *) &ibase;
  int n = sizeof(int);
  unsigned int hash = 0;
  for (i = 0; i < n; i++) {
    hash += str[i];
    hash += (hash << 10);
    hash ^= (hash >> 6);
  }
  str = (char *) coord;
  n = 3 * sizeof(double);
  for (i = 0; i < n; i++) {
    hash += str[i];
    hash += (hash << 10);
    hash ^= (hash >> 6);
  }
  hash += (hash << 3);
  hash ^= (hash >> 11);
  hash += (hash << 15);
  seed = hash & 0x7ffffff;
  if (!seed) seed = 1;
  for (i = 0; i < 5; i++) uniform();
  save = 0;
}
int RanPark::state()
{
  return seed;
}
