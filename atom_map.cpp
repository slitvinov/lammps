#include "atom.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include <cmath>
using namespace LAMMPS_NS;
#define EXTRA 1000
void Atom::map_delete() {
  memory->destroy(sametag);
  sametag = nullptr;
  max_same = 0;
  if (map_nhash) {
    delete[] map_bucket;
    delete[] map_hash;
    map_bucket = nullptr;
    map_hash = nullptr;
  }
  map_nhash = map_nbucket = 0;
}
