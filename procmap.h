#ifndef LMP_PROCMAP_H
#define LMP_PROCMAP_H
#include "pointers.h"
namespace LAMMPS_NS {
class ProcMap : protected Pointers {
public:
  ProcMap(class LAMMPS *);
  void onelevel_grid(int, int *, int *, int, int, int *, int *);
  void cart_map(int, int *, int *, int[3][2], int ***);
private:
  int procs_per_node;
  int procs_per_numa;
  int node_id;
  int nodegrid[3];
  int **cmap;
  int factor(int, int **);
  int combine_factors(int, int **, int, int **, int **);
  int cull_2d(int, int **, int);
  int cull_user(int, int **, int, int *);
  int cull_other(int, int **, int, int, int *, int *);
  int best_factors(int, int **, int *, int, int, int);
};
} // namespace LAMMPS_NS
#endif
