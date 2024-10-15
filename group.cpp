#include "group.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "exceptions.h"
#include "fix.h"
#include "force.h"
#include "input.h"
#include "math_eigen.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "region.h"
#include "tokenizer.h"
#include <cmath>
#include <cstring>
#include <map>
#include <utility>
using namespace LAMMPS_NS;
static constexpr int MAX_GROUP = 32;
static constexpr double EPSILON = 1.0e-6;
enum { NONE, TYPE, MOLECULE, ID };
enum { LT, LE, GT, GE, EQ, NEQ, BETWEEN };
#define BIG 1.0e20
Group::Group(LAMMPS *lmp) : Pointers(lmp) {
  MPI_Comm_rank(world, &me);
  names = new char *[MAX_GROUP];
  bitmask = new int[MAX_GROUP];
  inversemask = new int[MAX_GROUP];
  dynamic = new int[MAX_GROUP];
  for (int i = 0; i < MAX_GROUP; i++)
    names[i] = nullptr;
  for (int i = 0; i < MAX_GROUP; i++)
    bitmask[i] = 1 << i;
  for (int i = 0; i < MAX_GROUP; i++)
    inversemask[i] = bitmask[i] ^ ~0;
  for (int i = 0; i < MAX_GROUP; i++)
    dynamic[i] = 0;
  names[0] = utils::strdup("all");
  ngroup = 1;
}
Group::~Group() {
  for (int i = 0; i < MAX_GROUP; i++)
    delete[] names[i];
  delete[] names;
  delete[] bitmask;
  delete[] inversemask;
  delete[] dynamic;
}
int Group::find(const std::string &name) {
  for (int igroup = 0; igroup < MAX_GROUP; igroup++)
    if (names[igroup] && (name == names[igroup]))
      return igroup;
  return -1;
}
