#ifdef FIX_CLASS
// clang-format off
FixStyle(abf,FixABF);
// clang-format on
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

  char *inpfile;    // file to read rigid body attributes from

  int rstyle;       // SINGLE,MOLECULE,GROUP
  int setupflag;    // 1 if body properties are setup, else 0
  int earlyflag;    // 1 if forces/torques computed at post_force()

  int dimension;    // # of dimensions
  int nbody;        // # of rigid bodies
  int nlinear;      // # of linear rigid bodies
  int *nrigid;      // # of atoms in each rigid body
  int *mol2body;    // convert mol-ID to rigid body index
  int *body2mol;    // convert rigid body index to mol-ID
  int maxmol;       // size of mol2body = max mol-ID

  int *body;            // which body each atom is part of (-1 if none)
  double **displace;    // displacement of each atom in body coords

  double *masstotal;    // total mass of each rigid body
  double **xcm;         // coords of center-of-mass of each rigid body
  double **vcm;         // velocity of center-of-mass of each
  double **fcm;         // force on center-of-mass of each
  double **inertia;     // 3 principal components of inertia of each
  double **ex_space, **ey_space, **ez_space;
  // principal axes of each in space coords
  double **angmom;        // angular momentum of each in space coords
  double **omega;         // angular velocity of each in space coords
  double **torque;        // torque on each rigid body in space coords
  double **quat;          // quaternion of each rigid body
  imageint *imagebody;    // image flags of xcm of each rigid body

  double magn_mu_body[3];    // magnetic moment in reference of the ABF
  double **magn_mu;          // magnetic moment of the ABFs in laboratory frame
  double magn_omega;         // angular frequency of rotation of the magnetic field
  double magn_B_length;      // magnetic field magnitude
  double magn_B[3];          // magnetic field

  double **sum, **all;    // work vectors for each rigid body
  int **remapflag;        // PBC remap flags for each rigid body

  int reinitflag;    // 1 if re-initialize rigid bodies between runs

  imageint *xcmimage;    // internal image flags for atoms in rigid bodies
                         // set relative to in-box xcm of each body

  double tfactor;    // scale factor on temperature of rigid bodies

  int tstat_flag;    // NVT settings
  double t_start, t_stop, t_target;
  double t_period, t_freq;
  int t_chain, t_iter, t_order;

  int pstat_flag;    // NPT settings
  double p_start[3], p_stop[3];
  double p_period[3], p_freq[3];
  int p_flag[3];

  char *id_gravity;    // ID of fix gravity command to add gravity forces
  double *gvec;        // ptr to gravity vector inside the fix

  // map ABF position to direction perpendicular to magnetic field plane of
  // rotation
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

}    // namespace LAMMPS_NS

#endif
#endif
