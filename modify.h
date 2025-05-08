#ifndef LMP_MODIFY_H
#define LMP_MODIFY_H
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
  int nfix_restart_global;
  int nfix_restart_peratom;
  int nfix, maxfix;
  Fix **fix;
  int *fmask;
  Modify(class LAMMPS *);
  virtual void init();
  virtual void setup(int);
  virtual void initial_integrate(int);
  virtual void final_integrate();
  virtual void post_run();
  Fix *add_fix(int, char **, int trysuffix = 1);
  const std::vector<Fix *> &get_fix_list();

protected:
  int n_post_force, n_post_force_group;
  int *list_initial_integrate, *list_post_integrate;
  int *list_pre_exchange, *list_pre_neighbor, *list_post_neighbor;
  int *list_pre_force, *list_pre_reverse;
  int *list_post_force, *list_post_force_group;
  int *list_final_integrate, *list_end_of_step;
  int *list_energy_couple, *list_energy_global, *list_energy_atom;
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
