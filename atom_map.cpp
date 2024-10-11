#include "atom.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include <cmath>
using namespace LAMMPS_NS;
#define EXTRA 1000
void Atom::map_init(int check)
{
  int recreate = 0;
  if (check) recreate = map_style_set();
  if (map_style == MAP_ARRAY && map_tag_max > map_maxarray) recreate = 1;
  else if (map_style == MAP_HASH && nlocal+nghost > map_nhash) recreate = 1;
  if (!recreate) {
    if (map_style == MAP_ARRAY) {
      for (int i = 0; i <= map_tag_max; i++) map_array[i] = -1;
    } else {
      for (int i = 0; i < map_nbucket; i++) map_bucket[i] = -1;
      map_nused = 0;
      map_free = 0;
      for (int i = 0; i < map_nhash; i++) map_hash[i].next = i+1;
      if (map_nhash > 0) map_hash[map_nhash-1].next = -1;
    }
  } else {
    map_delete();
    if (map_style == MAP_ARRAY) {
      map_maxarray = map_tag_max;
      memory->create(map_array,map_maxarray+1,"atom:map_array");
      for (int i = 0; i <= map_tag_max; i++) map_array[i] = -1;
    } else {
      int nper = static_cast<int> (natoms/comm->nprocs);
      map_nhash = MAX(nper,nmax);
      map_nhash *= 2;
      map_nhash = MAX(map_nhash,1000);
      map_nbucket = next_prime(map_nhash);
      map_bucket = new int[map_nbucket];
      for (int i = 0; i < map_nbucket; i++) map_bucket[i] = -1;
      map_hash = new HashElem[map_nhash];
      map_nused = 0;
      map_free = 0;
      for (int i = 0; i < map_nhash; i++) map_hash[i].next = i+1;
      map_hash[map_nhash-1].next = -1;
    }
  }
}
void Atom::map_clear()
{
  if (map_style == MAP_ARRAY) {
    int nall = nlocal + nghost;
    for (int i = 0; i < nall; i++) {
      if (sametag) sametag[i] = -1;
      map_array[tag[i]] = -1;
    }
  } else {
    int previous,ibucket,index;
    tagint global;
    int nall = nlocal + nghost;
    for (int i = 0; i < nall; i++) {
      if (sametag) sametag[i] = -1;
      previous = -1;
      global = tag[i];
      ibucket = global % map_nbucket;
      index = map_bucket[ibucket];
      while (index > -1) {
        if (map_hash[index].global == global) break;
        previous = index;
        index = map_hash[index].next;
      }
      if (index == -1) continue;
      if (previous == -1) map_bucket[ibucket] = map_hash[index].next;
      else map_hash[previous].next = map_hash[index].next;
      map_hash[index].next = map_free;
      map_free = index;
      map_nused--;
    }
  }
}
void Atom::map_set()
{
  int nall = nlocal + nghost;
  if (map_style == MAP_ARRAY) {
    if (nall > max_same) {
      max_same = nall + EXTRA;
      memory->destroy(sametag);
      memory->create(sametag,max_same,"atom:sametag");
    }
    for (int i = nall-1; i >= 0 ; i--) {
      sametag[i] = map_array[tag[i]];
      map_array[tag[i]] = i;
    }
  } else {
    if (nall > map_nhash) map_init(0);
    if (nall > max_same) {
      max_same = nall + EXTRA;
      memory->destroy(sametag);
      memory->create(sametag,max_same,"atom:sametag");
    }
    int previous,ibucket,index;
    tagint global;
    for (int i = nall-1; i >= 0 ; i--) {
      sametag[i] = map_find_hash(tag[i]);
      previous = -1;
      global = tag[i];
      ibucket = global % map_nbucket;
      index = map_bucket[ibucket];
      while (index > -1) {
        if (map_hash[index].global == global) break;
        previous = index;
        index = map_hash[index].next;
      }
      if (index > -1) {
        map_hash[index].local = i;
        continue;
      }
      index = map_free;
      map_free = map_hash[map_free].next;
      if (previous == -1) map_bucket[ibucket] = index;
      else map_hash[previous].next = index;
      map_hash[index].global = global;
      map_hash[index].local = i;
      map_hash[index].next = -1;
      map_nused++;
    }
  }
}
void Atom::map_one(tagint global, int local)
{
  if (map_style == MAP_ARRAY) map_array[global] = local;
  else {
    int previous = -1;
    int ibucket = global % map_nbucket;
    int index = map_bucket[ibucket];
    while (index > -1) {
      if (map_hash[index].global == global) break;
      previous = index;
      index = map_hash[index].next;
    }
    if (index > -1) {
      map_hash[index].local = local;
      return;
    }
    index = map_free;
    map_free = map_hash[map_free].next;
    if (previous == -1) map_bucket[ibucket] = index;
    else map_hash[previous].next = index;
    map_hash[index].global = global;
    map_hash[index].local = local;
    map_hash[index].next = -1;
    map_nused++;
  }
}
int Atom::map_style_set()
{
  if (tag_enable == 0)
    error->all(FLERR,"Cannot create an atom map unless atoms have IDs");
  tagint max = -1;
  for (int i = 0; i < nlocal; i++) max = MAX(max,tag[i]);
  MPI_Allreduce(&max,&map_tag_max,1,MPI_LMP_TAGINT,MPI_MAX,world);
  int map_style_old = map_style;
  if (map_user == MAP_ARRAY || map_user == MAP_HASH) {
    map_style = map_user;
  } else {
    if (map_tag_max > 1000000) map_style = MAP_HASH;
    else map_style = MAP_ARRAY;
  }
  int recreate = 0;
  if (map_style != map_style_old) recreate = 1;
  return recreate;
}
void Atom::map_delete()
{
  memory->destroy(sametag);
  sametag = nullptr;
  max_same = 0;
  if (map_style == MAP_ARRAY) {
    memory->destroy(map_array);
    map_array = nullptr;
  } else {
    if (map_nhash) {
      delete [] map_bucket;
      delete [] map_hash;
      map_bucket = nullptr;
      map_hash = nullptr;
    }
    map_nhash = map_nbucket = 0;
  }
}
int Atom::map_find_hash(tagint global)
{
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
int Atom::next_prime(int n)
{
  int factor;
  int nprime = n+1;
  if (nprime % 2 == 0) nprime++;
  int root = static_cast<int> (sqrt(1.0*n)) + 2;
  while (nprime < MAXSMALLINT) {
    for (factor = 3; factor < root; factor++)
      if (nprime % factor == 0) break;
    if (factor == root) return nprime;
    nprime += 2;
  }
  return MAXSMALLINT;
}
