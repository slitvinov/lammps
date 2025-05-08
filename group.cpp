#include <map>
#include <set>
#include <unordered_set>
#include <cmath>
#include <cstring>
#include <map>
#include <utility>
#include <string>
#include <vector>
#include <mpi.h>
#include "lmptype.h"
#include "utils.h"
#include "lammps.h"
#include "pointers.h"
#include "lmptype.h"
#include "group.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "fix.h"
#include "force.h"
#include "memory.h"
using namespace LAMMPS_NS;
static constexpr int MAX_GROUP = 32;
static constexpr double EPSILON = 1.0e-6;
enum { NONE, TYPE, ID };
enum { LT, LE, GT, GE, EQ, NEQ, BETWEEN };
#define BIG 1.0e20
Group::Group(LAMMPS *lmp) : Pointers(lmp) {
  MPI_Comm_rank(world, &me);
  names = new char *[MAX_GROUP];
  bitmask = new int[MAX_GROUP];
  dynamic = new int[MAX_GROUP];
  for (int i = 0; i < MAX_GROUP; i++)
    names[i] = nullptr;
  for (int i = 0; i < MAX_GROUP; i++)
    bitmask[i] = 1 << i;
  for (int i = 0; i < MAX_GROUP; i++)
    dynamic[i] = 0;
  names[0] = utils::strdup("all");
  ngroup = 1;
}
int Group::find(const std::string &name) {
  for (int igroup = 0; igroup < MAX_GROUP; igroup++)
    if (names[igroup] && (name == names[igroup]))
      return igroup;
  return -1;
}
