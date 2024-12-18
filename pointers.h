#ifndef LMP_POINTERS_H
#define LMP_POINTERS_H
#include "fmt/format.h"
#include "lammps.h"
#include "lmptype.h"
#include "platform.h"
#include "utils.h"
#include <cstddef>
#include <cstdio>
#include <mpi.h>
#include <string>
#include <vector>
namespace LAMMPS_NS {
#define FLERR __FILE__, __LINE__
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#define MAX(A, B) ((A) > (B) ? (A) : (B))
enum ExecutionSpace { Host, Device };
template <class T> class MyPoolChunk;
template <class T> class MyPage;
class Pointers {
public:
  Pointers(LAMMPS *ptr)
      : lmp(ptr), memory(ptr->memory), error(ptr->error),
        universe(ptr->universe), input(ptr->input), atom(ptr->atom),
        update(ptr->update), neighbor(ptr->neighbor), comm(ptr->comm),
        domain(ptr->domain), force(ptr->force), modify(ptr->modify),
        group(ptr->group), world(ptr->world), infile(ptr->infile),
        screen(ptr->screen), logfile(ptr->logfile) {}
  virtual ~Pointers() = default;
  Pointers() = delete;
  Pointers(const Pointers &) = default;
  Pointers(Pointers &&) = delete;
  Pointers &operator=(const Pointers &) = delete;
  Pointers &operator=(Pointers &&) = delete;

protected:
  LAMMPS *lmp;
  Memory *&memory;
  Error *&error;
  Universe *&universe;
  Input *&input;
  Atom *&atom;
  Update *&update;
  Neighbor *&neighbor;
  Comm *&comm;
  Domain *&domain;
  Force *&force;
  Modify *&modify;
  Group *&group;
  MPI_Comm &world;
  FILE *&infile;
  FILE *&screen;
  FILE *&logfile;
};
} // namespace LAMMPS_NS
#endif
