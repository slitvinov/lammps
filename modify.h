#ifndef LMP_MODIFY_H
#define LMP_MODIFY_H
#include "pointers.h"
#include <map>
namespace LAMMPS_NS {
class Fix;
class Modify : protected Pointers {
  friend class Info;
  friend class FixSRP;
  friend class Respa;
  friend class RespaOMP;

public:
  int n_initial_integrate, n_post_integrate, n_pre_exchange;
  int n_pre_neighbor, n_post_neighbor;
  int n_pre_force, n_pre_reverse, n_post_force_any;
  int n_final_integrate, n_end_of_step;
  int n_energy_couple, n_energy_global, n_energy_atom;
  int n_initial_integrate_respa, n_post_integrate_respa;
  int n_pre_force_respa, n_post_force_respa_any, n_final_integrate_respa;
  int n_min_pre_exchange, n_min_pre_neighbor, n_min_post_neighbor;
  int n_min_pre_force, n_min_pre_reverse, n_min_post_force, n_min_energy;
  int restart_pbc_any;
  int nfix_restart_global;
  int nfix_restart_peratom;
  int nfix, maxfix;
  Fix **fix;
  int *fmask;
  Modify(class LAMMPS *);
  ~Modify() override;
  virtual void init();
  virtual void setup(int);
  virtual void setup_pre_exchange();
  virtual void setup_pre_neighbor();
  virtual void setup_post_neighbor();
  virtual void setup_pre_force(int);
  virtual void setup_pre_reverse(int, int);
  virtual void initial_integrate(int);
  virtual void post_integrate();
  virtual void pre_exchange();
  virtual void pre_neighbor();
  virtual void post_neighbor();
  virtual void pre_force(int);
  virtual void pre_reverse(int, int);
  virtual void post_force(int);
  virtual void final_integrate();
  virtual void end_of_step();
  virtual double energy_couple();
  virtual double energy_global();
  virtual void energy_atom(int, double *);
  virtual void post_run();
  virtual void setup_pre_force_respa(int, int);
  virtual void initial_integrate_respa(int, int, int);
  virtual void post_integrate_respa(int, int);
  virtual void pre_force_respa(int, int, int);
  virtual void post_force_respa(int, int, int);
  virtual void final_integrate_respa(int, int);
  virtual void min_pre_exchange();
  virtual void min_pre_neighbor();
  virtual void min_post_neighbor();
  virtual void min_pre_force(int);
  virtual void min_pre_reverse(int, int);
  virtual void min_post_force(int);
  virtual double min_energy(double *);
  virtual void min_store();
  virtual void min_step(double, double *);
  virtual void min_clearstore();
  virtual void min_pushstore();
  virtual void min_popstore();
  virtual double max_alpha(double *);
  virtual int min_dof();
  virtual int min_reset_ref();
  void reset_grid();
  Fix *add_fix(int, char **, int trysuffix = 1);
  Fix *add_fix(const std::string &, int trysuffix = 1);
  Fix *replace_fix(const char *, int, char **, int trysuffix = 1);
  Fix *replace_fix(const std::string &, const std::string &, int trysuffix = 1);
  void modify_fix(int, char **);
  void delete_fix(const std::string &);
  void delete_fix(int);
  int find_fix(const std::string &);
  Fix *get_fix_by_id(const std::string &) const;
  Fix *get_fix_by_index(int idx) const {
    return ((idx >= 0) && (idx < nfix)) ? fix[idx] : nullptr;
  }
  const std::vector<Fix *> get_fix_by_style(const std::string &) const;
  const std::vector<Fix *> &get_fix_list();
  int check_package(const char *);
  int check_rigid_group_overlap(int);
  int check_rigid_region_overlap(int, class Region *);
  int check_rigid_list_overlap(int *);

protected:
  int n_post_force, n_post_force_group, n_post_force_respa;
  int *list_initial_integrate, *list_post_integrate;
  int *list_pre_exchange, *list_pre_neighbor, *list_post_neighbor;
  int *list_pre_force, *list_pre_reverse;
  int *list_post_force, *list_post_force_group;
  int *list_final_integrate, *list_end_of_step;
  int *list_energy_couple, *list_energy_global, *list_energy_atom;
  int *list_initial_integrate_respa, *list_post_integrate_respa;
  int *list_pre_force_respa, *list_post_force_respa;
  int *list_final_integrate_respa;
  int *list_min_pre_exchange, *list_min_pre_neighbor, *list_min_post_neighbor;
  int *list_min_pre_force, *list_min_pre_reverse, *list_min_post_force;
  int *list_min_energy;
  int *end_of_step_every;
  int n_timeflag;
  int *list_timeflag;
  char **id_restart_global;
  char **style_restart_global;
  char **state_restart_global;
  int *used_restart_global;
  char **id_restart_peratom;
  char **style_restart_peratom;
  int *index_restart_peratom;
  int *used_restart_peratom;
  int index_permanent;
  std::vector<Fix *> fix_list;
  void list_init(int, int &, int *&);
  void list_init_end_of_step(int, int &, int *&);
  void list_init_energy_couple(int &, int *&);
  void list_init_energy_global(int &, int *&);
  void list_init_energy_atom(int &, int *&);
  void list_init_post_force_group(int &, int *&);
  void list_init_post_force_respa_group(int &, int *&);
  void list_init_dofflag(int &, int *&);

public:
  typedef Fix *(*FixCreator)(LAMMPS *, int, char **);
  typedef std::map<std::string, FixCreator> FixCreatorMap;
  FixCreatorMap *fix_map;

protected:
  void create_factories();
};
} // namespace LAMMPS_NS
#endif
