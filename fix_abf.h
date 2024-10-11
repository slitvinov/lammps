#ifdef FIX_CLASS
FixStyle(abf,FixABF);
#else
#ifndef LMP_FIX_ABF_H
#define LMP_FIX_ABF_H 
#include "fix.h"
namespace LAMMPS_NS {
class NeuralNet;
class FixABF : public Fix {
 public:
  FixABF(class LAMMPS *, int, char **);
  ~FixABF() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void initial_integrate(int) override;
  void post_force(int) override;
  void final_integrate() override;
  void write_restart_file(const char *) override;
  double compute_scalar() override;
  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void set_arrays(int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  void setup_pre_neighbor() override;
  void pre_neighbor() override;
  void deform(int) override;
  void enforce2d();
  void reset_dt() override;
  void zero_momentum() override;
  void zero_rotation() override;
  int modify_param(int, char **) override;
  void *extract(const char *, int &) override;
  double extract_ke();
  double extract_erotational();
  double compute_array(int, int) override;
 protected:
  int me, nprocs;
  double dtv, dtf, dtq;
  int triclinic;
  char *inpfile;
  int rstyle;
  int setupflag;
  int earlyflag;
  int dimension;
  int nbody;
  int nlinear;
  int *nrigid;
  int *mol2body;
  int *body2mol;
  int maxmol;
  int *body;
  double **displace;
  double *masstotal;
  double **xcm;
  double **vcm;
  double **fcm;
  double **inertia;
  double **ex_space, **ey_space, **ez_space;
  double **angmom;
  double **omega;
  double **torque;
  double **quat;
  imageint *imagebody;
  double magn_mu_body[3];
  double **magn_mu;
  double magn_omega;
  double magn_B_length;
  double magn_B[3];
  double **sum, **all;
  int **remapflag;
  int reinitflag;
  imageint *xcmimage;
  double tfactor;
  int tstat_flag;
  double t_start, t_stop, t_target;
  double t_period, t_freq;
  int t_chain, t_iter, t_order;
  int pstat_flag;
  double p_start[3], p_stop[3];
  double p_period[3], p_freq[3];
  int p_flag[3];
  char *id_gravity;
  double *gvec;
  NeuralNet *nn;
  MPI_Comm abfComm;
  bool isRankCloseToAbf;
  double abfCommMargin;
  int abfCommUpdatePeriod;
  double current_time();
  void image_shift();
  void set_xv();
  void set_v();
  void set_mu();
  void compute_B();
  void setup_bodies_static();
  void setup_bodies_dynamic();
  void compute_forces_and_torques();
  void readfile(int, double *, double **, double **, double **, imageint *, int *);
  void broadcastState();
  void updateAbfComm();
};
}
#endif
#endif
