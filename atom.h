#ifndef LMP_ATOM_H
#define LMP_ATOM_H
#include "pointers.h"
#include <map>
#include <set>
namespace LAMMPS_NS {
class AtomVec;
class Atom : protected Pointers {
public:
  char *atom_style;
  AtomVec *avec;
  enum { DOUBLE, INT, BIGINT };
  enum { GROW = 0, RESTART = 1, BORDER = 2 };
  enum { ATOMIC = 0 };
  enum { ATOM = 0 };
  enum { NUMERIC = 0, LABELS = 1 };
  enum { MAP_NONE = 0, MAP_ARRAY = 1, MAP_HASH = 2, MAP_YES = 3 };
  bigint natoms;
  int nlocal, nghost;
  int nmax;
  int tag_enable;
  int ntypes;
  int firstgroup;
  int nfirst;
  char *firstgroupname;
  tagint *tag;
  int *type, *mask;
  imageint *image;
  double **x, **v, **f;
  double *rmass;
  double *q, **mu;
  double *radius;
  double **omega, **angmom, **torque;
  double **quat;
  double *temperature, *heatflow;
  double *vfrac, *s0;
  double **x0;
  double **sp, **fm, **fm_long;
  int *spin;
  double *eradius, *ervel, *erforce;
  double *ervelforce;
  double **cs, **csforce, **vforce;
  int *etag;
  tagint *id5p;
  double **cc, **cc_flux;
  double *edpd_temp, *edpd_flux;
  double *vest_temp;
  double *edpd_cv;
  int cc_species;
  double *contact_radius;
  double **smd_data_9;
  double **smd_stress;
  double *eff_plastic_strain;
  double *eff_plastic_strain_rate;
  double *damage;
  double *rho, *drho, *esph, *desph, *cv;
  double **vest;
  int types_style;
  int peri_flag, electron_flag;
  int wavepacket_flag, sph_flag;
  int q_flag, mu_flag;
  int rmass_flag, radius_flag, omega_flag, torque_flag, angmom_flag, quat_flag;
  int temperature_flag, heatflow_flag;
  int vfrac_flag, spin_flag, eradius_flag, ervel_flag, erforce_flag;
  int cs_flag, csforce_flag, vforce_flag, ervelforce_flag, etag_flag;
  int rho_flag, esph_flag, cv_flag, vest_flag;
  int dpd_flag, edpd_flag, tdpd_flag;
  int mesont_flag;
  struct PerAtom {
    std::string name;
    void *address;
    void *address_length;
    int *address_maxcols;
    int datatype;
    int cols;
    int collength;
    int threadflag;
  };
  std::vector<PerAtom> peratom;
  int **ivector, ***iarray;
  double **dvector, ***darray;
  int *icols, *dcols;
  char **ivname, **dvname, **ianame, **daname;
  int nivector, ndvector, niarray, ndarray;
  double **extra;
  double *mass;
  int *mass_setflag;
  int nextra_grow, nextra_restart, nextra_border;
  int *extra_grow, *extra_restart, *extra_border;
  int nextra_grow_max, nextra_restart_max;
  int nextra_border_max;
  int nextra_store;
  int map_style;
  int map_user;
  tagint map_tag_max;
  std::set<tagint> *unique_tags;
  int sortfreq;
  bigint nextsort;
  double userbinsize;
  int *sametag;
  bool reset_image_flag[3];
  typedef AtomVec *(*AtomVecCreator)(LAMMPS *);
  typedef std::map<std::string, AtomVecCreator> AtomVecCreatorMap;
  AtomVecCreatorMap *avec_map;
  Atom(class LAMMPS *);
  ~Atom() override;
  void peratom_create();
  void add_peratom(const std::string &, void *, int, int, int threadflag = 0);
  void create_avec(const std::string &, int, char **, int);
  virtual AtomVec *new_avec(const std::string &);
  void init();
  void setup();
  AtomVec *style_match(const char *);
  void modify_params(int, char **);
  void tag_check();
  void tag_extend();
  int tag_consecutive();
  int parse_data(const char *);
  virtual void allocate_type_arrays();
  void set_mass(const char *, int, const char *, int, int, int *);
  void set_mass(const char *, int, int, double);
  void set_mass(const char *, int, int, char **);
  void set_mass(double *);
  void check_mass(const char *, int);
  void first_reorder();
  virtual void sort();
  void add_callback(int);
  void delete_callback(const char *, int);
  void update_callback(int);
  int find_custom(const char *, int &, int &);
  virtual int add_custom(const char *, int, int);
  virtual void sync_modify(ExecutionSpace, unsigned int, unsigned int) {}
  inline int *get_map_array() { return map_array; };
  inline int get_map_size() { return map_tag_max + 1; };
  inline int get_max_same() { return max_same; };
  inline int get_map_maxarray() { return map_maxarray + 1; };
  int memcheck(const char *) { return 1; }
  inline int map(tagint global) {
    if (map_style == 1)
      return map_array[global];
    else if (map_style == 2)
      return map_find_hash(global);
    else
      return -1;
  };
  virtual void map_init(int check = 1);
  virtual void map_clear();
  virtual void map_set();
  void map_one(tagint, int);
  int map_style_set();
  virtual void map_delete();
  int map_find_hash(tagint);

protected:
  int *map_array;
  int map_maxarray;
  struct HashElem {
    tagint global;
    int local;
    int next;
  };
  int map_nhash;
  int map_nused;
  int map_free;
  int map_nbucket;
  int *map_bucket;
  HashElem *map_hash;
  int max_same;
  int nbins;
  int nbinx, nbiny, nbinz;
  int maxbin;
  int maxnext;
  int *binhead;
  int *next;
  int *permute;
  double bininvx, bininvy, bininvz;
  double bboxlo[3], bboxhi[3];
  void set_atomflag_defaults();
  void setup_sort_bins();
  int next_prime(int);
};
} // namespace LAMMPS_NS
#endif
