#ifndef LMP_NPAIR_H
#define LMP_NPAIR_H
namespace LAMMPS_NS {
class NPair : protected Pointers {
public:
  int istyle;
  class NBin *nb;
  bigint last_build;
  double cutoff_custom;
  NPair(class LAMMPS *);
  void post_constructor(class NeighRequest *);
  virtual void copy_neighbor_info();
  void build_setup();
  virtual void build(class NeighList *) = 0;

protected:
  double **mycutneighsq;
  int includegroup;
  int exclude;
  double skin;
  double **cutneighsq;
  double **cutneighghostsq;
  double cut_inner_sq;
  double cut_middle_sq;
  double cut_middle_inside_sq;
  double *bboxlo, *bboxhi;
  int ncollections;
  double **cutcollectionsq;
  int nex_type;
  int *ex1_type, *ex2_type;
  int **ex_type;
  int nex_group;
  int *ex1_group, *ex2_group;
  int *ex1_bit, *ex2_bit;
  int nex_mol;
  int *ex_mol_bit;
  int *ex_mol_group;
  int *ex_mol_intra;
  int *special_flag;
  int nbinx, nbiny, nbinz;
  int mbins;
  int mbinx, mbiny, mbinz;
  int mbinxlo, mbinylo, mbinzlo;
  double bininvx, bininvy, bininvz;
  int *atom2bin, *bins;
  int *binhead;
  int *nbinx_multi, *nbiny_multi, *nbinz_multi;
  int *mbins_multi;
  int *mbinx_multi, *mbiny_multi, *mbinz_multi;
  int *mbinxlo_multi, *mbinylo_multi, *mbinzlo_multi;
  double *bininvx_multi, *bininvy_multi, *bininvz_multi;
  double **distsq_multi_old;
  virtual void copy_bin_info();
  inline int find_special(const tagint *list, const int *nspecial,
                          const tagint tag) const {
    const int n1 = nspecial[0];
    const int n2 = nspecial[1];
    const int n3 = nspecial[2];
    for (int i = 0; i < n3; i++) {
      if (list[i] == tag) {
        if (i < n1) {
          if (special_flag[1] == 0)
            return -1;
          else if (special_flag[1] == 1)
            return 0;
          else
            return 1;
        } else if (i < n2) {
          if (special_flag[2] == 0)
            return -1;
          else if (special_flag[2] == 1)
            return 0;
          else
            return 2;
        } else {
          if (special_flag[3] == 0)
            return -1;
          else if (special_flag[3] == 1)
            return 0;
          else
            return 3;
        }
      }
    }
    return 0;
  };
  int copymode;
  ExecutionSpace execution_space;
};
} // namespace LAMMPS_NS
#endif
