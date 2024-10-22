#include "atom.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include <cmath>
using namespace LAMMPS_NS;
#define EXTRA 1000
void Atom::map_clear() {
  if (map_style == MAP_ARRAY) {
    int nall = nlocal + nghost;
    for (int i = 0; i < nall; i++) {
      if (sametag)
        sametag[i] = -1;
      map_array[tag[i]] = -1;
    }
  } else {
    int previous, ibucket, index;
    tagint global;
    int nall = nlocal + nghost;
    for (int i = 0; i < nall; i++) {
      if (sametag)
        sametag[i] = -1;
      previous = -1;
      global = tag[i];
      ibucket = global % map_nbucket;
      index = map_bucket[ibucket];
      while (index > -1) {
        if (map_hash[index].global == global)
          break;
        previous = index;
        index = map_hash[index].next;
      }
      if (index == -1)
        continue;
      if (previous == -1)
        map_bucket[ibucket] = map_hash[index].next;
      else
        map_hash[previous].next = map_hash[index].next;
      map_hash[index].next = map_free;
      map_free = index;
      map_nused--;
    }
  }
}
void Atom::map_delete() {
  memory->destroy(sametag);
  sametag = nullptr;
  max_same = 0;
  if (map_style == MAP_ARRAY) {
    memory->destroy(map_array);
    map_array = nullptr;
  } else {
    if (map_nhash) {
      delete[] map_bucket;
      delete[] map_hash;
      map_bucket = nullptr;
      map_hash = nullptr;
    }
    map_nhash = map_nbucket = 0;
  }
}
int Atom::map_find_hash(tagint global) {
  int local = -1;
  int index = map_bucket[global % map_nbucket];
  while (index > -1) {
    if (map_hash[index].global == global) {
      local = map_hash[index].local;
      break;
    }
    index = map_hash[index].next;
  }
  return local;
}
