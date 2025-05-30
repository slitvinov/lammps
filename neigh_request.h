#ifndef LMP_NEIGH_REQUEST_H
#define LMP_NEIGH_REQUEST_H
namespace LAMMPS_NS {
class NeighRequest : protected Pointers {
  friend class Neighbor;
  friend class NBin;
  friend class NeighList;
  friend class NPair;

protected:
  void *requestor;
  int requestor_instance;
  int id;
  int pair;
  int fix;
  int compute;
  int command;
  int neigh;
  int half;
  int full;
  int occasional;
  int newton;
  int ghost;
  int size;
  int history;
  int granonesided;
  int respainner;
  int respamiddle;
  int respaouter;
  int bond;
  int omp;
  int ssa;
  int cut;
  double cutoff;
  int skip;
  int *iskip;
  int **ijskip;
  const char *command_style;
  int skiplist;
  int off2on;
  int copy;
  int trim;
  int copylist;
  int halffull;
  int halffulllist;
  int unique;
  int index_bin;
  int index_pair;

public:
  NeighRequest(class LAMMPS *);
  NeighRequest(class LAMMPS *, void *, int);
  NeighRequest(NeighRequest *);
  ~NeighRequest() override;
  void copy_request(NeighRequest *, int);
  int get_size() const { return size; }
  void *get_requestor() const { return requestor; }
};
} // namespace LAMMPS_NS
#endif
