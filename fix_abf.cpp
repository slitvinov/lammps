#include "fix_abf.h"

#include "abf_const.h"
#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "math_const.h"
#include "math_eigen.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "tokenizer.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;
using namespace AbfConst;

namespace LAMMPS_NS {

struct LinearLayer {
  int nx, ny;
  std::vector<double> A, b;
};

std::vector<LinearLayer> readLayers(const char *path)
{
  FILE *f = fopen(path, "r");

  if (f == nullptr) {
    printf("Could not open %s\n", path);
    exit(1);
  }

  std::vector<LinearLayer> layers;

  while (true) {
    LinearLayer l;

    if (2 != fscanf(f, "linear %d %d\n", &l.nx, &l.ny)) break;
    l.A.resize(l.nx * l.ny);
    l.b.resize(l.nx);

    fscanf(f, "A\n");
    for (int i = 0; i < l.nx * l.ny; ++i) fscanf(f, "%lg ", &l.A[i]);

    fscanf(f, "b\n");
    for (int i = 0; i < l.nx; ++i) fscanf(f, "%lg ", &l.b[i]);

    layers.push_back(std::move(l));
  }
  fclose(f);

  return layers;
}

class NeuralNet {
 public:
  NeuralNet(const char *path) : NeuralNet(readLayers(path)) {}

  NeuralNet(std::vector<LinearLayer> layers) : _layers(std::move(layers))
  {
    int maxw = 0;
    for (size_t i = 0; i < _layers.size(); ++i) {
      maxw = std::max(maxw, _layers[i].nx);
      maxw = std::max(maxw, _layers[i].ny);
    }

    _work0.resize(maxw);
    _work1.resize(maxw);
  }

  ~NeuralNet() = default;

  void eval(const double *in, double *out) const
  {
    const double *x = in;
    const size_t K = _layers.size();

    for (size_t k = 0; k < K; ++k) {
      double *y;
      if (k + 1 == K)
        y = out;
      else
        y = _work0.data();

      const LinearLayer &l = _layers[k];

      for (size_t i = 0; i < l.nx; ++i) {
        y[i] = l.b[i];

        for (size_t j = 0; j < l.ny; ++j) y[i] += l.A[i * l.ny + j] * x[j];
      }

      if (k + 1 != K) {
        for (size_t i = 0; i < l.nx; ++i) y[i] = std::tanh(y[i]);

        x = y;
        std::swap(_work0, _work1);
      }
    }
  }

 private:
  std::vector<LinearLayer> _layers;
  mutable std::vector<double> _work0, _work1;
};

}    // namespace LAMMPS_NS

/* ---------------------------------------------------------------------- */

FixABF::FixABF(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), inpfile(nullptr), nrigid(nullptr), mol2body(nullptr), body2mol(nullptr),
    body(nullptr), displace(nullptr), masstotal(nullptr), xcm(nullptr), vcm(nullptr), fcm(nullptr),
    inertia(nullptr), ex_space(nullptr), ey_space(nullptr), ez_space(nullptr), angmom(nullptr),
    omega(nullptr), torque(nullptr), quat(nullptr), imagebody(nullptr), magn_mu_body{0.0, 0.0, 0.0},
    magn_mu(nullptr), magn_omega(0.0), magn_B_length(0.0), magn_B{0.0, 0.0, 0.0}, sum(nullptr),
    all(nullptr), remapflag(nullptr), xcmimage(nullptr), id_gravity(nullptr), nn(nullptr),
    abfComm(MPI_COMM_NULL), isRankCloseToAbf(true), abfCommMargin(0.0), abfCommUpdatePeriod(0)
{
  int i, ibody;

  scalar_flag = 1;
  extscalar = 0;
  time_integrate = 1;
  rigid_flag = 1;
  virial_global_flag = virial_peratom_flag = 1;
  thermo_virial = 1;
  create_attribute = 1;
  dof_flag = 1;
  centroidstressflag = CENTROID_NOTAVAIL;

  MPI_Comm_rank(world, &me);
  MPI_Comm_size(world, &nprocs);

  // perform initial allocation of atom-based arrays
  // register with Atom class

  body = nullptr;
  xcmimage = nullptr;
  displace = nullptr;
  FixABF::grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);

  // parse args for rigid body specification
  // set nbody and body[i] for each atom

  if (narg < 4) error->all(FLERR, "Illegal fix abf command");
  int iarg;

  mol2body = nullptr;
  body2mol = nullptr;

  // single abf
  // nbody = 1
  // all atoms in fix group are part of body

  if (strcmp(arg[3], "single") == 0) {
    rstyle = SINGLE;
    iarg = 4;
    nbody = 1;

    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++) {
      body[i] = -1;
      if (mask[i] & groupbit) body[i] = 0;
    }

  } else {
    error->all(FLERR, "Illegal fix abf command: only 'single' is supported");
  }

  // error check on nbody

  if (nbody == 0) error->all(FLERR, "No rigid bodies defined");

  // create all nbody-length arrays

  memory->create(nrigid, nbody, "abf:nrigid");
  memory->create(masstotal, nbody, "abf:masstotal");
  memory->create(xcm, nbody, 3, "abf:xcm");
  memory->create(vcm, nbody, 3, "abf:vcm");
  memory->create(fcm, nbody, 3, "abf:fcm");
  memory->create(inertia, nbody, 3, "abf:inertia");
  memory->create(ex_space, nbody, 3, "abf:ex_space");
  memory->create(ey_space, nbody, 3, "abf:ey_space");
  memory->create(ez_space, nbody, 3, "abf:ez_space");
  memory->create(angmom, nbody, 3, "abf:angmom");
  memory->create(omega, nbody, 3, "abf:omega");
  memory->create(torque, nbody, 3, "abf:torque");
  memory->create(quat, nbody, 4, "abf:quat");
  memory->create(imagebody, nbody, "abf:imagebody");
  memory->create(magn_mu, nbody, 3, "abf:mu");

  memory->create(sum, nbody, 6, "abf:sum");
  memory->create(all, nbody, 6, "abf:all");
  memory->create(remapflag, nbody, 4, "abf:remapflag");

  array_flag = 1;
  size_array_rows = nbody;
  size_array_cols = 15;
  global_freq = 1;
  extarray = 0;

  // number of linear rigid bodies is counted later

  nlinear = 0;

  // parse optional args

  reinitflag = 1;

  tstat_flag = 0;
  pstat_flag = 0;
  t_chain = 10;
  t_iter = 1;
  t_order = 3;

  inpfile = nullptr;
  id_gravity = nullptr;

  dimension = domain->dimension;

  for (i = 0; i < 3; i++) {
    p_start[i] = p_stop[i] = p_period[i] = 0.0;
    p_flag[i] = 0;
  }

  while (iarg < narg) {
    if (strcmp(arg[iarg], "mu") == 0) {
      if (iarg + 4 > narg) error->all(FLERR, "Illegal fix abf command: mu");

      magn_mu_body[0] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      magn_mu_body[1] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      magn_mu_body[2] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);

      iarg += 4;

    } else if (strcmp(arg[iarg], "omega") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix abf command: omega");

      magn_omega = utils::numeric(FLERR, arg[iarg + 1], false, lmp);

      iarg += 2;

    } else if (strcmp(arg[iarg], "B") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix abf command: B");

      magn_B_length = utils::numeric(FLERR, arg[iarg + 1], false, lmp);

      iarg += 2;

    } else if (strcmp(arg[iarg], "NN") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix abf command: NN");

      nn = new NeuralNet(arg[iarg + 1]);

      iarg += 2;

    } else if (strcmp(arg[iarg], "abf_comm") == 0) {
      if (iarg + 3 > narg) error->all(FLERR, "Illegal fix abf command: abf_comm");

      abfCommMargin = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      abfCommUpdatePeriod = utils::inumeric(FLERR, arg[iarg + 2], false, lmp);

      iarg += 3;

    } else if (strcmp(arg[iarg], "temp") == 0) {
      if (iarg + 4 > narg) error->all(FLERR, "Illegal fix rigid command");
      if (!utils::strmatch(style, "^rigid/n.t")) error->all(FLERR, "Illegal fix rigid command");
      tstat_flag = 1;
      t_start = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      t_stop = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      t_period = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      iarg += 4;

    } else if (strcmp(arg[iarg], "iso") == 0) {
      if (iarg + 4 > narg) error->all(FLERR, "Illegal fix rigid command");
      if (!utils::strmatch(style, "^rigid/np.")) error->all(FLERR, "Illegal fix rigid command");
      p_start[0] = p_start[1] = p_start[2] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      p_stop[0] = p_stop[1] = p_stop[2] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      p_period[0] = p_period[1] = p_period[2] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      p_flag[0] = p_flag[1] = p_flag[2] = 1;
      if (dimension == 2) {
        p_start[2] = p_stop[2] = p_period[2] = 0.0;
        p_flag[2] = 0;
      }
      iarg += 4;

    } else if (strcmp(arg[iarg], "aniso") == 0) {
      if (iarg + 4 > narg) error->all(FLERR, "Illegal fix rigid command");
      if (!utils::strmatch(style, "^rigid/np.")) error->all(FLERR, "Illegal fix rigid command");
      p_start[0] = p_start[1] = p_start[2] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      p_stop[0] = p_stop[1] = p_stop[2] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      p_period[0] = p_period[1] = p_period[2] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      p_flag[0] = p_flag[1] = p_flag[2] = 1;
      if (dimension == 2) {
        p_start[2] = p_stop[2] = p_period[2] = 0.0;
        p_flag[2] = 0;
      }
      iarg += 4;

    } else if (strcmp(arg[iarg], "x") == 0) {
      if (iarg + 4 > narg) error->all(FLERR, "Illegal fix rigid command");
      if (!utils::strmatch(style, "^rigid/np.")) error->all(FLERR, "Illegal fix rigid command");
      p_start[0] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      p_stop[0] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      p_period[0] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      p_flag[0] = 1;
      iarg += 4;

    } else if (strcmp(arg[iarg], "y") == 0) {
      if (iarg + 4 > narg) error->all(FLERR, "Illegal fix rigid command");
      if (!utils::strmatch(style, "^rigid/np.")) error->all(FLERR, "Illegal fix rigid command");
      p_start[1] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      p_stop[1] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      p_period[1] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      p_flag[1] = 1;
      iarg += 4;

    } else if (strcmp(arg[iarg], "z") == 0) {
      if (iarg + 4 > narg) error->all(FLERR, "Illegal fix rigid command");
      if (!utils::strmatch(style, "^rigid/np.")) error->all(FLERR, "Illegal fix rigid command");
      p_start[2] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      p_stop[2] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      p_period[2] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      p_flag[2] = 1;
      iarg += 4;

    } else if (strcmp(arg[iarg], "tparam") == 0) {
      if (iarg + 4 > narg) error->all(FLERR, "Illegal fix rigid command");
      if (!utils::strmatch(style, "^rigid/n.t")) error->all(FLERR, "Illegal fix rigid command");
      t_chain = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      t_iter = utils::inumeric(FLERR, arg[iarg + 2], false, lmp);
      t_order = utils::inumeric(FLERR, arg[iarg + 3], false, lmp);
      iarg += 4;

    } else if (strcmp(arg[iarg], "infile") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix rigid command");
      delete[] inpfile;
      inpfile = utils::strdup(arg[iarg + 1]);
      restart_file = 1;
      reinitflag = 0;
      iarg += 2;

    } else if (strcmp(arg[iarg], "reinit") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix rigid command");
      reinitflag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;

    } else if (strcmp(arg[iarg], "gravity") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix rigid command");
      delete[] id_gravity;
      id_gravity = utils::strdup(arg[iarg + 1]);
      iarg += 2;

    } else
      error->all(FLERR, "Illegal fix rigid command");
  }

  // set pstat_flag

  pstat_flag = 0;
  for (i = 0; i < 3; i++)
    if (p_flag[i]) pstat_flag = 1;

  // initialize vector output quantities in case accessed before run

  for (i = 0; i < nbody; i++) {
    xcm[i][0] = xcm[i][1] = xcm[i][2] = 0.0;
    vcm[i][0] = vcm[i][1] = vcm[i][2] = 0.0;
    fcm[i][0] = fcm[i][1] = fcm[i][2] = 0.0;
    torque[i][0] = torque[i][1] = torque[i][2] = 0.0;
  }

  // nrigid[n] = # of atoms in Nth rigid body
  // error if one or zero atoms

  int *ncount = new int[nbody];
  for (ibody = 0; ibody < nbody; ibody++) ncount[ibody] = 0;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++)
    if (body[i] >= 0) ncount[body[i]]++;

  MPI_Allreduce(ncount, nrigid, nbody, MPI_INT, MPI_SUM, world);
  delete[] ncount;

  for (ibody = 0; ibody < nbody; ibody++)
    if (nrigid[ibody] <= 1) error->all(FLERR, "One or zero atoms in rigid body");

  // wait to setup bodies until first init() using current atom properties

  setupflag = 0;

  // compute per body forces and torques at final_integrate() by default

  earlyflag = 0;

  // print statistics

  int nsum = 0;
  for (ibody = 0; ibody < nbody; ibody++) nsum += nrigid[ibody];

  if (me == 0) utils::logmesg(lmp, "  {} rigid bodies with {} atoms\n", nbody, nsum);
}

/* ---------------------------------------------------------------------- */

FixABF::~FixABF()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id, Atom::GROW);

  delete[] inpfile;
  delete[] id_gravity;

  memory->destroy(mol2body);
  memory->destroy(body2mol);

  // delete locally stored per-atom arrays

  memory->destroy(body);
  memory->destroy(xcmimage);
  memory->destroy(displace);

  // delete nbody-length arrays

  memory->destroy(nrigid);
  memory->destroy(masstotal);
  memory->destroy(xcm);
  memory->destroy(vcm);
  memory->destroy(fcm);
  memory->destroy(inertia);
  memory->destroy(ex_space);
  memory->destroy(ey_space);
  memory->destroy(ez_space);
  memory->destroy(angmom);
  memory->destroy(omega);
  memory->destroy(torque);
  memory->destroy(quat);
  memory->destroy(imagebody);
  memory->destroy(magn_mu);

  memory->destroy(sum);
  memory->destroy(all);
  memory->destroy(remapflag);

  if (nn) delete nn;

  if (abfComm != MPI_COMM_NULL) MPI_Comm_free(&abfComm);
}

/* ---------------------------------------------------------------------- */

int FixABF::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= POST_FORCE;    // magnetic moment
  mask |= PRE_NEIGHBOR;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixABF::init()
{
  triclinic = domain->triclinic;

  // warn if more than one rigid fix
  // if earlyflag, warn if any post-force fixes come after a rigid fix

  int count = 0;
  for (auto &ifix : modify->get_fix_list())
    if (ifix->rigid_flag) count++;
  if (count > 1 && me == 0) error->warning(FLERR, "More than one fix rigid");

  if (earlyflag) {
    bool rflag = false;
    for (auto &ifix : modify->get_fix_list()) {
      if (ifix->rigid_flag) rflag = true;
      if ((comm->me == 0) && rflag && (ifix->setmask() & POST_FORCE) && !ifix->rigid_flag)
        error->warning(FLERR, "Fix {} with ID {} alters forces after fix rigid", ifix->style,
                       ifix->id);
    }
  }

  // warn if body properties are read from inpfile
  //   and the gravity keyword is not set and a gravity fix exists
  // this could mean body particles are overlapped
  //   and gravity is not applied correctly

  if (inpfile && !id_gravity) {
    if (modify->get_fix_by_style("^gravity").size() > 0)
      if (comm->me == 0)
        error->warning(FLERR,
                       "Gravity may not be correctly applied to rigid "
                       "bodies if they consist of overlapped particles");
  }

  //  error if a fix changing the box comes before rigid fix

  bool boxflag = false;
  for (auto &ifix : modify->get_fix_list()) {
    if (boxflag && utils::strmatch(ifix->style, "^rigid"))
      error->all(FLERR, "Rigid fixes must come before any box changing fix");
    if (ifix->box_change) boxflag = true;
  }

  // add gravity forces based on gravity vector from fix

  if (id_gravity) {
    auto ifix = modify->get_fix_by_id(id_gravity);
    if (!ifix) error->all(FLERR, "Fix rigid cannot find fix gravity ID {}", id_gravity);
    if (!utils::strmatch(ifix->style, "^gravity"))
      error->all(FLERR, "Fix rigid gravity fix ID {} is not a gravity fix style", id_gravity);
    int tmp;
    gvec = (double *) ifix->extract("gvec", tmp);
  }

  // timestep info

  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  dtq = 0.5 * update->dt;

  // setup rigid bodies, using current atom info. if reinitflag is not set,
  // do the initialization only once, b/c properties may not be re-computable
  // especially if overlapping particles.
  //   do not do dynamic init if read body properties from inpfile.
  // this is b/c the inpfile defines the static and dynamic properties and may
  // not be computable if contain overlapping particles.
  //   setup_bodies_static() reads inpfile itself

  if (reinitflag || !setupflag) {
    setup_bodies_static();
    if (!inpfile) setup_bodies_dynamic();
    setupflag = 1;
  }

  // temperature scale factor

  double ndof = 6 * nbody;
  ndof -= nlinear;
  if (ndof > 0.0)
    tfactor = force->mvv2e / (ndof * force->boltz);
  else
    tfactor = 0.0;

  // abf communicator
  if (abfCommUpdatePeriod > 0) { updateAbfComm(); }
}

/* ----------------------------------------------------------------------
   invoke pre_neighbor() to ensure body xcmimage flags are reset
     needed if Verlet::setup::pbc() has remapped/migrated atoms for 2nd run
------------------------------------------------------------------------- */

void FixABF::setup_pre_neighbor()
{
  pre_neighbor();
}

/* ----------------------------------------------------------------------
   compute initial fcm and torque on bodies, also initial virial
   reset all particle velocities to be consistent with vcm and omega
------------------------------------------------------------------------- */

void FixABF::setup(int vflag)
{
  int i, n, ibody;

  // fcm = force on center-of-mass of each rigid body

  double **f = atom->f;
  int nlocal = atom->nlocal;

  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;

  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];
    sum[ibody][0] += f[i][0];
    sum[ibody][1] += f[i][1];
    sum[ibody][2] += f[i][2];
  }

  MPI_Allreduce(sum[0], all[0], 6 * nbody, MPI_DOUBLE, MPI_SUM, world);

  for (ibody = 0; ibody < nbody; ibody++) {
    fcm[ibody][0] = all[ibody][0];
    fcm[ibody][1] = all[ibody][1];
    fcm[ibody][2] = all[ibody][2];
  }

  // torque = torque on each rigid body

  double **x = atom->x;

  double dx, dy, dz;
  double unwrap[3];

  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;

  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    domain->unmap(x[i], xcmimage[i], unwrap);
    dx = unwrap[0] - xcm[ibody][0];
    dy = unwrap[1] - xcm[ibody][1];
    dz = unwrap[2] - xcm[ibody][2];

    sum[ibody][0] += dy * f[i][2] - dz * f[i][1];
    sum[ibody][1] += dz * f[i][0] - dx * f[i][2];
    sum[ibody][2] += dx * f[i][1] - dy * f[i][0];
  }

  MPI_Allreduce(sum[0], all[0], 6 * nbody, MPI_DOUBLE, MPI_SUM, world);

  for (ibody = 0; ibody < nbody; ibody++) {
    torque[ibody][0] = all[ibody][0];
    torque[ibody][1] = all[ibody][1];
    torque[ibody][2] = all[ibody][2];
  }

  // virial setup before call to set_v

  v_init(vflag);

  // set velocities from angmom & omega

  for (ibody = 0; ibody < nbody; ibody++)
    MathExtra::angmom_to_omega(angmom[ibody], ex_space[ibody], ey_space[ibody], ez_space[ibody],
                               inertia[ibody], omega[ibody]);

  set_v();

  // guesstimate virial as 2x the set_v contribution

  if (vflag_global)
    for (n = 0; n < 6; n++) virial[n] *= 2.0;
  if (vflag_atom) {
    for (i = 0; i < nlocal; i++)
      for (n = 0; n < 6; n++) vatom[i][n] *= 2.0;
  }
}

/* ---------------------------------------------------------------------- */

void FixABF::initial_integrate(int vflag)
{
  double dtfm;

  for (int ibody = 0; ibody < nbody; ibody++) {

    // update vcm by 1/2 step

    dtfm = dtf / masstotal[ibody];
    vcm[ibody][0] += dtfm * fcm[ibody][0];
    vcm[ibody][1] += dtfm * fcm[ibody][1];
    vcm[ibody][2] += dtfm * fcm[ibody][2];

    // update xcm by full step

    xcm[ibody][0] += dtv * vcm[ibody][0];
    xcm[ibody][1] += dtv * vcm[ibody][1];
    xcm[ibody][2] += dtv * vcm[ibody][2];

    // update angular momentum by 1/2 step

    angmom[ibody][0] += dtf * torque[ibody][0];
    angmom[ibody][1] += dtf * torque[ibody][1];
    angmom[ibody][2] += dtf * torque[ibody][2];

    // compute omega at 1/2 step from angmom at 1/2 step and current q
    // update quaternion a full step via Richardson iteration
    // returns new normalized quaternion, also updated omega at 1/2 step
    // update ex,ey,ez to reflect new quaternion

    MathExtra::angmom_to_omega(angmom[ibody], ex_space[ibody], ey_space[ibody], ez_space[ibody],
                               inertia[ibody], omega[ibody]);
    MathExtra::richardson(quat[ibody], angmom[ibody], omega[ibody], inertia[ibody], dtq);
    MathExtra::q_to_exyz(quat[ibody], ex_space[ibody], ey_space[ibody], ez_space[ibody]);
  }

  // virial setup before call to set_xv

  v_init(vflag);

  // set coords/orient and velocity/rotation of atoms in rigid bodies
  // from quarternion and omega

  set_xv();
}

/* ----------------------------------------------------------------------
   called from FixEnforce2d post_force() for 2d problems
   zero all body values that should be zero for 2d model
------------------------------------------------------------------------- */

void FixABF::enforce2d()
{
  for (int ibody = 0; ibody < nbody; ibody++) {
    xcm[ibody][2] = 0.0;
    vcm[ibody][2] = 0.0;
    fcm[ibody][2] = 0.0;
    torque[ibody][0] = 0.0;
    torque[ibody][1] = 0.0;
    angmom[ibody][0] = 0.0;
    angmom[ibody][1] = 0.0;
    omega[ibody][0] = 0.0;
    omega[ibody][1] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

void FixABF::compute_forces_and_torques()
{

  if (abfCommUpdatePeriod > 0 && update->ntimestep % abfCommUpdatePeriod == 0) {
    broadcastState();
    updateAbfComm();
  }

  int i, ibody;

  // sum over atoms to get force and torque on rigid body

  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;

  double dx, dy, dz;
  double unwrap[3];

  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;

  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    sum[ibody][0] += f[i][0];
    sum[ibody][1] += f[i][1];
    sum[ibody][2] += f[i][2];

    domain->unmap(x[i], xcmimage[i], unwrap);
    dx = unwrap[0] - xcm[ibody][0];
    dy = unwrap[1] - xcm[ibody][1];
    dz = unwrap[2] - xcm[ibody][2];

    sum[ibody][3] += dy * f[i][2] - dz * f[i][1];
    sum[ibody][4] += dz * f[i][0] - dx * f[i][2];
    sum[ibody][5] += dx * f[i][1] - dy * f[i][0];
  }

  if (abfCommUpdatePeriod > 0) {
    if (isRankCloseToAbf) {
      MPI_Allreduce(sum[0], all[0], 6 * nbody, MPI_DOUBLE, MPI_SUM, abfComm);
    } else {
      for (ibody = 0; ibody < nbody; ++ibody) {
        for (i = 0; i < 6; ++i) {
          if (sum[ibody][i] != 0)
            error->one(FLERR, "Non-zero data was not communicated. check abfComm update");
        }
      }
    }

  } else {
    MPI_Allreduce(sum[0], all[0], 6 * nbody, MPI_DOUBLE, MPI_SUM, world);
  }

  for (ibody = 0; ibody < nbody; ibody++) {
    fcm[ibody][0] = all[ibody][0];
    fcm[ibody][1] = all[ibody][1];
    fcm[ibody][2] = all[ibody][2];
    torque[ibody][0] = all[ibody][3];
    torque[ibody][1] = all[ibody][4];
    torque[ibody][2] = all[ibody][5];
  }

  // add magnetic torque
  set_mu();
  compute_B();
  double T[3];
  for (ibody = 0; ibody < nbody; ibody++) {
    MathExtra::cross3(magn_mu[ibody], magn_B, T);
    torque[ibody][0] += T[0];
    torque[ibody][1] += T[1];
    torque[ibody][2] += T[2];
  }

  // add gravity force to COM of each body

  if (id_gravity) {
    for (ibody = 0; ibody < nbody; ibody++) {
      fcm[ibody][0] += gvec[0] * masstotal[ibody];
      fcm[ibody][1] += gvec[1] * masstotal[ibody];
      fcm[ibody][2] += gvec[2] * masstotal[ibody];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixABF::post_force(int /*vflag*/)
{
  if (earlyflag) compute_forces_and_torques();
}

/* ---------------------------------------------------------------------- */

void FixABF::final_integrate()
{
  int ibody;
  double dtfm;

  if (!earlyflag) compute_forces_and_torques();

  // update vcm and angmom

  for (ibody = 0; ibody < nbody; ibody++) {

    // update vcm by 1/2 step

    dtfm = dtf / masstotal[ibody];
    vcm[ibody][0] += dtfm * fcm[ibody][0];
    vcm[ibody][1] += dtfm * fcm[ibody][1];
    vcm[ibody][2] += dtfm * fcm[ibody][2];

    // update angular momentum by 1/2 step

    angmom[ibody][0] += dtf * torque[ibody][0];
    angmom[ibody][1] += dtf * torque[ibody][1];
    angmom[ibody][2] += dtf * torque[ibody][2];

    MathExtra::angmom_to_omega(angmom[ibody], ex_space[ibody], ey_space[ibody], ez_space[ibody],
                               inertia[ibody], omega[ibody]);
  }

  // set velocity/rotation of atoms in rigid bodies
  // virial is already setup from initial_integrate

  set_v();
}

/* ----------------------------------------------------------------------
   remap xcm of each rigid body back into periodic simulation box
   done during pre_neighbor so will be after call to pbc()
     and after fix_deform::pre_exchange() may have flipped box
   use domain->remap() in case xcm is far away from box
     due to first-time definition of rigid body in setup_bodies_static()
     or due to box flip
   also adjust imagebody = rigid body image flags, due to xcm remap
   also reset body xcmimage flags of all atoms in bodies
   xcmimage flags are relative to xcm so that body can be unwrapped
   if don't do this, would need xcm to move with true image flags
     then a body could end up very far away from box
     set_xv() will then compute huge displacements every step to
       reset coords of all body atoms to be back inside the box,
       ditto for triclinic box flip, which causes numeric problems
------------------------------------------------------------------------- */

void FixABF::pre_neighbor()
{
  for (int ibody = 0; ibody < nbody; ibody++) domain->remap(xcm[ibody], imagebody[ibody]);
  image_shift();
}

double FixABF::current_time()
{
  return update->atime + (update->ntimestep - update->atimestep) * update->dt;
}

/* ----------------------------------------------------------------------
   reset body xcmimage flags of atoms in bodies
   xcmimage flags are relative to xcm so that body can be unwrapped
   xcmimage = true image flag - imagebody flag
------------------------------------------------------------------------- */

void FixABF::image_shift()
{
  int ibody;
  imageint tdim, bdim, xdim[3];

  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    tdim = image[i] & IMGMASK;
    bdim = imagebody[ibody] & IMGMASK;
    xdim[0] = IMGMAX + tdim - bdim;
    tdim = (image[i] >> IMGBITS) & IMGMASK;
    bdim = (imagebody[ibody] >> IMGBITS) & IMGMASK;
    xdim[1] = IMGMAX + tdim - bdim;
    tdim = image[i] >> IMG2BITS;
    bdim = imagebody[ibody] >> IMG2BITS;
    xdim[2] = IMGMAX + tdim - bdim;

    xcmimage[i] = (xdim[2] << IMG2BITS) | (xdim[1] << IMGBITS) | xdim[0];
  }
}

/* ----------------------------------------------------------------------
   adjust xcm of each rigid body due to box deformation
   called by various fixes that change box size/shape
   flag = 0/1 means map from box to lamda coords or vice versa
------------------------------------------------------------------------- */

void FixABF::deform(int flag)
{
  if (flag == 0)
    for (int ibody = 0; ibody < nbody; ibody++) domain->x2lamda(xcm[ibody], xcm[ibody]);
  else
    for (int ibody = 0; ibody < nbody; ibody++) domain->lamda2x(xcm[ibody], xcm[ibody]);
}

/* ----------------------------------------------------------------------
   set space-frame coords and velocity of each atom in each rigid body
   set orientation and rotation of extended particles
   x = Q displace + Xcm, mapped back to periodic box
   v = Vcm + (W cross (x - Xcm))
------------------------------------------------------------------------- */

void FixABF::set_xv()
{
  int ibody;
  int xbox, ybox, zbox;
  double x0, x1, x2, v0, v1, v2, fc0, fc1, fc2, massone;
  double xy, xz, yz;
  double ione[3], exone[3], eyone[3], ezone[3], vr[6], p[3][3];

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  if (triclinic) {
    xy = domain->xy;
    xz = domain->xz;
    yz = domain->yz;
  }

  // set x and v of each atom

  for (int i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    xbox = (xcmimage[i] & IMGMASK) - IMGMAX;
    ybox = (xcmimage[i] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (xcmimage[i] >> IMG2BITS) - IMGMAX;

    // save old positions and velocities for virial

    if (evflag) {
      if (triclinic == 0) {
        x0 = x[i][0] + xbox * xprd;
        x1 = x[i][1] + ybox * yprd;
        x2 = x[i][2] + zbox * zprd;
      } else {
        x0 = x[i][0] + xbox * xprd + ybox * xy + zbox * xz;
        x1 = x[i][1] + ybox * yprd + zbox * yz;
        x2 = x[i][2] + zbox * zprd;
      }
      v0 = v[i][0];
      v1 = v[i][1];
      v2 = v[i][2];
    }

    // x = displacement from center-of-mass, based on body orientation
    // v = vcm + omega around center-of-mass

    MathExtra::matvec(ex_space[ibody], ey_space[ibody], ez_space[ibody], displace[i], x[i]);

    v[i][0] = omega[ibody][1] * x[i][2] - omega[ibody][2] * x[i][1] + vcm[ibody][0];
    v[i][1] = omega[ibody][2] * x[i][0] - omega[ibody][0] * x[i][2] + vcm[ibody][1];
    v[i][2] = omega[ibody][0] * x[i][1] - omega[ibody][1] * x[i][0] + vcm[ibody][2];

    // add center of mass to displacement
    // map back into periodic box via xbox,ybox,zbox
    // for triclinic, add in box tilt factors as well

    if (triclinic == 0) {
      x[i][0] += xcm[ibody][0] - xbox * xprd;
      x[i][1] += xcm[ibody][1] - ybox * yprd;
      x[i][2] += xcm[ibody][2] - zbox * zprd;
    } else {
      x[i][0] += xcm[ibody][0] - xbox * xprd - ybox * xy - zbox * xz;
      x[i][1] += xcm[ibody][1] - ybox * yprd - zbox * yz;
      x[i][2] += xcm[ibody][2] - zbox * zprd;
    }

    // virial = unwrapped coords dotted into body constraint force
    // body constraint force = implied force due to v change minus f external
    // assume f does not include forces internal to body
    // 1/2 factor b/c final_integrate contributes other half
    // assume per-atom contribution is due to constraint force on that atom

    if (evflag) {
      if (rmass)
        massone = rmass[i];
      else
        massone = mass[type[i]];
      fc0 = massone * (v[i][0] - v0) / dtf - f[i][0];
      fc1 = massone * (v[i][1] - v1) / dtf - f[i][1];
      fc2 = massone * (v[i][2] - v2) / dtf - f[i][2];

      vr[0] = 0.5 * x0 * fc0;
      vr[1] = 0.5 * x1 * fc1;
      vr[2] = 0.5 * x2 * fc2;
      vr[3] = 0.5 * x0 * fc1;
      vr[4] = 0.5 * x0 * fc2;
      vr[5] = 0.5 * x1 * fc2;

      v_tally(1, &i, 1.0, vr);
    }
  }
}

/* ----------------------------------------------------------------------
   set space-frame velocity of each atom in a rigid body
   set omega and angmom of extended particles
   v = Vcm + (W cross (x - Xcm))
------------------------------------------------------------------------- */

void FixABF::set_v()
{
  int xbox, ybox, zbox;
  double x0, x1, x2, v0, v1, v2, fc0, fc1, fc2, massone;
  double xy, xz, yz;
  double ione[3], exone[3], eyone[3], ezone[3], delta[3], vr[6];

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  if (triclinic) {
    xy = domain->xy;
    xz = domain->xz;
    yz = domain->yz;
  }

  // set v of each atom

  for (int i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    const int ibody = body[i];

    MathExtra::matvec(ex_space[ibody], ey_space[ibody], ez_space[ibody], displace[i], delta);

    // save old velocities for virial

    if (evflag) {
      v0 = v[i][0];
      v1 = v[i][1];
      v2 = v[i][2];
    }

    v[i][0] = omega[ibody][1] * delta[2] - omega[ibody][2] * delta[1] + vcm[ibody][0];
    v[i][1] = omega[ibody][2] * delta[0] - omega[ibody][0] * delta[2] + vcm[ibody][1];
    v[i][2] = omega[ibody][0] * delta[1] - omega[ibody][1] * delta[0] + vcm[ibody][2];

    // virial = unwrapped coords dotted into body constraint force
    // body constraint force = implied force due to v change minus f external
    // assume f does not include forces internal to body
    // 1/2 factor b/c initial_integrate contributes other half
    // assume per-atom contribution is due to constraint force on that atom

    if (evflag) {
      if (rmass)
        massone = rmass[i];
      else
        massone = mass[type[i]];
      fc0 = massone * (v[i][0] - v0) / dtf - f[i][0];
      fc1 = massone * (v[i][1] - v1) / dtf - f[i][1];
      fc2 = massone * (v[i][2] - v2) / dtf - f[i][2];

      xbox = (xcmimage[i] & IMGMASK) - IMGMAX;
      ybox = (xcmimage[i] >> IMGBITS & IMGMASK) - IMGMAX;
      zbox = (xcmimage[i] >> IMG2BITS) - IMGMAX;

      if (triclinic == 0) {
        x0 = x[i][0] + xbox * xprd;
        x1 = x[i][1] + ybox * yprd;
        x2 = x[i][2] + zbox * zprd;
      } else {
        x0 = x[i][0] + xbox * xprd + ybox * xy + zbox * xz;
        x1 = x[i][1] + ybox * yprd + zbox * yz;
        x2 = x[i][2] + zbox * zprd;
      }

      vr[0] = 0.5 * x0 * fc0;
      vr[1] = 0.5 * x1 * fc1;
      vr[2] = 0.5 * x2 * fc2;
      vr[3] = 0.5 * x0 * fc1;
      vr[4] = 0.5 * x0 * fc2;
      vr[5] = 0.5 * x1 * fc2;

      v_tally(1, &i, 1.0, vr);
    }
  }
}

void FixABF::set_mu()
{
  for (int i = 0; i < nbody; ++i) {
    MathExtra::matvec(ex_space[i], ey_space[i], ez_space[i], magn_mu_body, magn_mu[i]);
  }
}

void FixABF::compute_B()
{
  const double t = current_time();
  const double ex[3] = {1.0, 0.0, 0.0};
  double direction[3] = {1.0, 0.0, 0.0};

  if (nn) {
    if (nbody != 1) error->all(FLERR, "Illegal use of NN: assume a single ABF.");
    nn->eval(xcm[0], direction);
    MathExtra::norm3(direction);
  }

  double q[4], axis[3];

  const double angle = std::acos(MathExtra::dot3(ex, direction));

  MathExtra::cross3(ex, direction, axis);
  if (MathExtra::len3(axis) == 0.0)
    axis[1] = 1.0;
  else
    MathExtra::norm3(axis);

  MathExtra::axisangle_to_quat(axis, angle, q);

  double B_[3] = {0.0, magn_B_length * std::cos(magn_omega * t),
                  magn_B_length * std::sin(magn_omega * t)};

  MathExtra::quatrotvec(q, B_, magn_B);
}

/* ----------------------------------------------------------------------
   one-time initialization of static rigid body attributes
   sets extended flags, masstotal, center-of-mass
   sets Cartesian and diagonalized inertia tensor
   sets body image flags
   may read some properties from inpfile
------------------------------------------------------------------------- */

void FixABF::setup_bodies_static()
{
  int i, ibody;

  double **mu = atom->mu;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *ellipsoid = atom->ellipsoid;
  int *line = atom->line;
  int *tri = atom->tri;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  // set body xcmimage flags = true image flags

  imageint *image = atom->image;
  for (i = 0; i < nlocal; i++)
    if (body[i] >= 0)
      xcmimage[i] = image[i];
    else
      xcmimage[i] = 0;

  // compute masstotal & center-of-mass of each rigid body
  // error if image flag is not 0 in a non-periodic dim

  double **x = atom->x;

  int *periodicity = domain->periodicity;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double xy = domain->xy;
  double xz = domain->xz;
  double yz = domain->yz;

  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;
  int xbox, ybox, zbox;
  double massone, xunwrap, yunwrap, zunwrap;

  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    xbox = (xcmimage[i] & IMGMASK) - IMGMAX;
    ybox = (xcmimage[i] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (xcmimage[i] >> IMG2BITS) - IMGMAX;
    if (rmass)
      massone = rmass[i];
    else
      massone = mass[type[i]];

    if ((xbox && !periodicity[0]) || (ybox && !periodicity[1]) || (zbox && !periodicity[2]))
      error->one(FLERR,
                 "Fix rigid atom has non-zero image flag "
                 "in a non-periodic dimension");

    if (triclinic == 0) {
      xunwrap = x[i][0] + xbox * xprd;
      yunwrap = x[i][1] + ybox * yprd;
      zunwrap = x[i][2] + zbox * zprd;
    } else {
      xunwrap = x[i][0] + xbox * xprd + ybox * xy + zbox * xz;
      yunwrap = x[i][1] + ybox * yprd + zbox * yz;
      zunwrap = x[i][2] + zbox * zprd;
    }

    sum[ibody][0] += xunwrap * massone;
    sum[ibody][1] += yunwrap * massone;
    sum[ibody][2] += zunwrap * massone;
    sum[ibody][3] += massone;
  }

  MPI_Allreduce(sum[0], all[0], 6 * nbody, MPI_DOUBLE, MPI_SUM, world);

  for (ibody = 0; ibody < nbody; ibody++) {
    masstotal[ibody] = all[ibody][3];
    xcm[ibody][0] = all[ibody][0] / masstotal[ibody];
    xcm[ibody][1] = all[ibody][1] / masstotal[ibody];
    xcm[ibody][2] = all[ibody][2] / masstotal[ibody];
  }

  // set vcm, angmom = 0.0 in case inpfile is used
  // and doesn't overwrite all body's values
  // since setup_bodies_dynamic() will not be called

  for (ibody = 0; ibody < nbody; ibody++) {
    vcm[ibody][0] = vcm[ibody][1] = vcm[ibody][2] = 0.0;
    angmom[ibody][0] = angmom[ibody][1] = angmom[ibody][2] = 0.0;
  }

  // set rigid body image flags to default values

  for (ibody = 0; ibody < nbody; ibody++)
    imagebody[ibody] = ((imageint) IMGMAX << IMG2BITS) | ((imageint) IMGMAX << IMGBITS) | IMGMAX;

  // overwrite masstotal, center-of-mass, image flags with file values
  // inbody[i] = 0/1 if Ith rigid body is initialized by file

  int *inbody;
  if (inpfile) {
    // must call it here so it doesn't override read in data but
    // initialize bodies whose dynamic settings not set in inpfile

    setup_bodies_dynamic();

    memory->create(inbody, nbody, "abf:inbody");
    for (ibody = 0; ibody < nbody; ibody++) inbody[ibody] = 0;
    readfile(0, masstotal, xcm, vcm, angmom, imagebody, inbody);
  }

  // remap the xcm of each body back into simulation box
  //   and reset body and atom xcmimage flags via pre_neighbor()

  pre_neighbor();

  // compute 6 moments of inertia of each body in Cartesian reference frame
  // dx,dy,dz = coords relative to center-of-mass
  // symmetric 3x3 inertia tensor stored in Voigt notation as 6-vector

  double dx, dy, dz;

  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;

  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    xbox = (xcmimage[i] & IMGMASK) - IMGMAX;
    ybox = (xcmimage[i] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (xcmimage[i] >> IMG2BITS) - IMGMAX;

    if (triclinic == 0) {
      xunwrap = x[i][0] + xbox * xprd;
      yunwrap = x[i][1] + ybox * yprd;
      zunwrap = x[i][2] + zbox * zprd;
    } else {
      xunwrap = x[i][0] + xbox * xprd + ybox * xy + zbox * xz;
      yunwrap = x[i][1] + ybox * yprd + zbox * yz;
      zunwrap = x[i][2] + zbox * zprd;
    }

    dx = xunwrap - xcm[ibody][0];
    dy = yunwrap - xcm[ibody][1];
    dz = zunwrap - xcm[ibody][2];

    if (rmass)
      massone = rmass[i];
    else
      massone = mass[type[i]];

    sum[ibody][0] += massone * (dy * dy + dz * dz);
    sum[ibody][1] += massone * (dx * dx + dz * dz);
    sum[ibody][2] += massone * (dx * dx + dy * dy);
    sum[ibody][3] -= massone * dy * dz;
    sum[ibody][4] -= massone * dx * dz;
    sum[ibody][5] -= massone * dx * dy;
  }

  MPI_Allreduce(sum[0], all[0], 6 * nbody, MPI_DOUBLE, MPI_SUM, world);

  // overwrite Cartesian inertia tensor with file values

  if (inpfile) readfile(1, nullptr, all, nullptr, nullptr, nullptr, inbody);

  // diagonalize inertia tensor for each body via Jacobi rotations
  // inertia = 3 eigenvalues = principal moments of inertia
  // evectors and exzy_space = 3 evectors = principal axes of rigid body

  int ierror;
  double cross[3];
  double tensor[3][3], evectors[3][3];

  for (ibody = 0; ibody < nbody; ibody++) {
    tensor[0][0] = all[ibody][0];
    tensor[1][1] = all[ibody][1];
    tensor[2][2] = all[ibody][2];
    tensor[1][2] = tensor[2][1] = all[ibody][3];
    tensor[0][2] = tensor[2][0] = all[ibody][4];
    tensor[0][1] = tensor[1][0] = all[ibody][5];

    ierror = MathEigen::jacobi3(tensor, inertia[ibody], evectors);
    if (ierror) error->all(FLERR, "Insufficient Jacobi rotations for rigid body");

    ex_space[ibody][0] = evectors[0][0];
    ex_space[ibody][1] = evectors[1][0];
    ex_space[ibody][2] = evectors[2][0];
    ey_space[ibody][0] = evectors[0][1];
    ey_space[ibody][1] = evectors[1][1];
    ey_space[ibody][2] = evectors[2][1];
    ez_space[ibody][0] = evectors[0][2];
    ez_space[ibody][1] = evectors[1][2];
    ez_space[ibody][2] = evectors[2][2];

    // if any principal moment < scaled EPSILON, set to 0.0

    double max;
    max = MAX(inertia[ibody][0], inertia[ibody][1]);
    max = MAX(max, inertia[ibody][2]);

    if (inertia[ibody][0] < EPSILON * max) inertia[ibody][0] = 0.0;
    if (inertia[ibody][1] < EPSILON * max) inertia[ibody][1] = 0.0;
    if (inertia[ibody][2] < EPSILON * max) inertia[ibody][2] = 0.0;

    // enforce 3 evectors as a right-handed coordinate system
    // flip 3rd vector if needed

    MathExtra::cross3(ex_space[ibody], ey_space[ibody], cross);
    if (MathExtra::dot3(cross, ez_space[ibody]) < 0.0) MathExtra::negate3(ez_space[ibody]);

    // create initial quaternion

    MathExtra::exyz_to_q(ex_space[ibody], ey_space[ibody], ez_space[ibody], quat[ibody]);
  }

  // displace = initial atom coords in basis of principal axes
  // set displace = 0.0 for atoms not in any rigid body
  // for extended particles, set their orientation wrt to rigid body

  double qc[4], delta[3];
  double *quatatom;
  double theta_body;

  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) {
      displace[i][0] = displace[i][1] = displace[i][2] = 0.0;
      continue;
    }

    ibody = body[i];

    xbox = (xcmimage[i] & IMGMASK) - IMGMAX;
    ybox = (xcmimage[i] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (xcmimage[i] >> IMG2BITS) - IMGMAX;

    if (triclinic == 0) {
      xunwrap = x[i][0] + xbox * xprd;
      yunwrap = x[i][1] + ybox * yprd;
      zunwrap = x[i][2] + zbox * zprd;
    } else {
      xunwrap = x[i][0] + xbox * xprd + ybox * xy + zbox * xz;
      yunwrap = x[i][1] + ybox * yprd + zbox * yz;
      zunwrap = x[i][2] + zbox * zprd;
    }

    delta[0] = xunwrap - xcm[ibody][0];
    delta[1] = yunwrap - xcm[ibody][1];
    delta[2] = zunwrap - xcm[ibody][2];
    MathExtra::transpose_matvec(ex_space[ibody], ey_space[ibody], ez_space[ibody], delta,
                                displace[i]);
  }

  // test for valid principal moments & axes
  // recompute moments of inertia around new axes
  // 3 diagonal moments should equal principal moments
  // 3 off-diagonal moments should be 0.0
  // extended particles may contribute extra terms to moments of inertia

  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;

  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];
    if (rmass)
      massone = rmass[i];
    else
      massone = mass[type[i]];

    sum[ibody][0] += massone * (displace[i][1] * displace[i][1] + displace[i][2] * displace[i][2]);
    sum[ibody][1] += massone * (displace[i][0] * displace[i][0] + displace[i][2] * displace[i][2]);
    sum[ibody][2] += massone * (displace[i][0] * displace[i][0] + displace[i][1] * displace[i][1]);
    sum[ibody][3] -= massone * displace[i][1] * displace[i][2];
    sum[ibody][4] -= massone * displace[i][0] * displace[i][2];
    sum[ibody][5] -= massone * displace[i][0] * displace[i][1];
  }

  MPI_Allreduce(sum[0], all[0], 6 * nbody, MPI_DOUBLE, MPI_SUM, world);

  // error check that re-computed moments of inertia match diagonalized ones
  // do not do test for bodies with params read from inpfile

  double norm;
  for (ibody = 0; ibody < nbody; ibody++) {
    if (inpfile && inbody[ibody]) continue;
    if (inertia[ibody][0] == 0.0) {
      if (fabs(all[ibody][0]) > TOLERANCE) error->all(FLERR, "Fix abf: Bad principal moments");
    } else {
      if (fabs((all[ibody][0] - inertia[ibody][0]) / inertia[ibody][0]) > TOLERANCE)
        error->all(FLERR, "Fix abf: Bad principal moments");
    }
    if (inertia[ibody][1] == 0.0) {
      if (fabs(all[ibody][1]) > TOLERANCE) error->all(FLERR, "Fix abf: Bad principal moments");
    } else {
      if (fabs((all[ibody][1] - inertia[ibody][1]) / inertia[ibody][1]) > TOLERANCE)
        error->all(FLERR, "Fix abf: Bad principal moments");
    }
    if (inertia[ibody][2] == 0.0) {
      if (fabs(all[ibody][2]) > TOLERANCE) error->all(FLERR, "Fix abf: Bad principal moments");
    } else {
      if (fabs((all[ibody][2] - inertia[ibody][2]) / inertia[ibody][2]) > TOLERANCE)
        error->all(FLERR, "Fix abf: Bad principal moments");
    }
    norm = (inertia[ibody][0] + inertia[ibody][1] + inertia[ibody][2]) / 3.0;
    if (fabs(all[ibody][3] / norm) > TOLERANCE || fabs(all[ibody][4] / norm) > TOLERANCE ||
        fabs(all[ibody][5] / norm) > TOLERANCE)
      error->all(FLERR, "Fix abf: Bad principal moments");
  }

  if (inpfile) memory->destroy(inbody);
}

/* ----------------------------------------------------------------------
   one-time initialization of dynamic rigid body attributes
   set vcm and angmom, computed explicitly from constituent particles
   not done if body properties read from file, e.g. for overlapping particles
------------------------------------------------------------------------- */

void FixABF::setup_bodies_dynamic()
{
  int i, ibody;
  double massone, radone;

  // vcm = velocity of center-of-mass of each rigid body
  // angmom = angular momentum of each rigid body

  double **x = atom->x;
  double **v = atom->v;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double dx, dy, dz;
  double unwrap[3];

  for (ibody = 0; ibody < nbody; ibody++)
    for (i = 0; i < 6; i++) sum[ibody][i] = 0.0;

  for (i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    ibody = body[i];

    if (rmass)
      massone = rmass[i];
    else
      massone = mass[type[i]];

    sum[ibody][0] += v[i][0] * massone;
    sum[ibody][1] += v[i][1] * massone;
    sum[ibody][2] += v[i][2] * massone;

    domain->unmap(x[i], xcmimage[i], unwrap);
    dx = unwrap[0] - xcm[ibody][0];
    dy = unwrap[1] - xcm[ibody][1];
    dz = unwrap[2] - xcm[ibody][2];

    sum[ibody][3] += dy * massone * v[i][2] - dz * massone * v[i][1];
    sum[ibody][4] += dz * massone * v[i][0] - dx * massone * v[i][2];
    sum[ibody][5] += dx * massone * v[i][1] - dy * massone * v[i][0];
  }

  MPI_Allreduce(sum[0], all[0], 6 * nbody, MPI_DOUBLE, MPI_SUM, world);

  // normalize velocity of COM

  for (ibody = 0; ibody < nbody; ibody++) {
    vcm[ibody][0] = all[ibody][0] / masstotal[ibody];
    vcm[ibody][1] = all[ibody][1] / masstotal[ibody];
    vcm[ibody][2] = all[ibody][2] / masstotal[ibody];
    angmom[ibody][0] = all[ibody][3];
    angmom[ibody][1] = all[ibody][4];
    angmom[ibody][2] = all[ibody][5];
  }
}

/* ----------------------------------------------------------------------
   read per rigid body info from user-provided file
   which = 0 to read everything except 6 moments of inertia
   which = 1 to read 6 moments of inertia
   flag inbody = 0 for bodies whose info is read from file
   nlines = # of lines of rigid body info
   one line = rigid-ID mass xcm ycm zcm ixx iyy izz ixy ixz iyz
              vxcm vycm vzcm lx ly lz ix iy iz
------------------------------------------------------------------------- */

void FixABF::readfile(int which, double *vec, double **array1, double **array2, double **array3,
                      imageint *ivec, int *inbody)
{
  int nchunk, id, eofflag, xbox, ybox, zbox;
  int nlines;
  FILE *fp;
  char *eof, *start, *next, *buf;
  char line[MAXLINE];

  // open file and read and parse first non-empty, non-comment line containing the number of bodies
  if (me == 0) {
    fp = fopen(inpfile, "r");
    if (fp == nullptr)
      error->one(FLERR, "Cannot open fix rigid infile {}: {}", inpfile, utils::getsyserror());
    while (true) {
      eof = fgets(line, MAXLINE, fp);
      if (eof == nullptr) error->one(FLERR, "Unexpected end of fix rigid infile");
      start = &line[strspn(line, " \t\n\v\f\r")];
      if (*start != '\0' && *start != '#') break;
    }
    nlines = utils::inumeric(FLERR, utils::trim(line), true, lmp);
    if (which == 0)
      utils::logmesg(lmp, "Reading rigid body data for {} bodies from file {}\n", nlines, inpfile);
    if (nlines == 0) fclose(fp);
  }
  MPI_Bcast(&nlines, 1, MPI_INT, 0, world);

  // empty file with 0 lines is needed to trigger initial restart file
  // generation when no infile was previously used.

  if (nlines == 0)
    return;
  else if (nlines < 0)
    error->all(FLERR, "Fix rigid infile has incorrect format");

  auto buffer = new char[CHUNK * MAXLINE];
  int nread = 0;
  while (nread < nlines) {
    nchunk = MIN(nlines - nread, CHUNK);
    eofflag = utils::read_lines_from_file(fp, nchunk, MAXLINE, buffer, me, world);
    if (eofflag) error->all(FLERR, "Unexpected end of fix rigid infile");

    buf = buffer;
    next = strchr(buf, '\n');
    *next = '\0';
    int nwords = utils::count_words(utils::trim_comment(buf));
    *next = '\n';

    if (nwords != ATTRIBUTE_PERBODY)
      error->all(FLERR, "Incorrect rigid body format in fix rigid file");

    // loop over lines of rigid body attributes
    // tokenize the line into values
    // id = rigid body ID
    // use ID as-is for SINGLE, as mol-ID for MOLECULE, as-is for GROUP
    // for which = 0, store all but inertia in vecs and arrays
    // for which = 1, store inertia tensor array, invert 3,4,5 values to Voigt

    for (int i = 0; i < nchunk; i++) {
      next = strchr(buf, '\n');
      *next = '\0';

      try {
        ValueTokenizer values(buf);
        id = values.next_int();
        if (rstyle == MOLECULE) {
          if (id <= 0 || id > maxmol)
            throw TokenizerException("invalid rigid molecule ID ", std::to_string(id));
          id = mol2body[id];
        } else
          id--;

        if (id < 0 || id >= nbody)
          throw TokenizerException("invalid_rigid body ID ", std::to_string(id + 1));

        inbody[id] = 1;

        if (which == 0) {
          vec[id] = values.next_double();
          array1[id][0] = values.next_double();
          array1[id][1] = values.next_double();
          array1[id][2] = values.next_double();
          values.skip(6);
          array2[id][0] = values.next_double();
          array2[id][1] = values.next_double();
          array2[id][2] = values.next_double();
          array3[id][0] = values.next_double();
          array3[id][1] = values.next_double();
          array3[id][2] = values.next_double();
          xbox = values.next_int();
          ybox = values.next_int();
          zbox = values.next_int();
          ivec[id] = ((imageint)(xbox + IMGMAX) & IMGMASK) |
              (((imageint)(ybox + IMGMAX) & IMGMASK) << IMGBITS) |
              (((imageint)(zbox + IMGMAX) & IMGMASK) << IMG2BITS);
        } else {
          values.skip(4);
          array1[id][0] = values.next_double();
          array1[id][1] = values.next_double();
          array1[id][2] = values.next_double();
          array1[id][5] = values.next_double();
          array1[id][4] = values.next_double();
          array1[id][3] = values.next_double();
        }
      } catch (TokenizerException &e) {
        error->all(FLERR, "Invalid fix rigid infile: {}", e.what());
      }
      buf = next + 1;
    }
    nread += nchunk;
  }

  if (me == 0) fclose(fp);
  delete[] buffer;
}

/* ----------------------------------------------------------------------
   write out restart info for mass, COM, inertia tensor, image flags to file
   identical format to inpfile option, so info can be read in when restarting
   only proc 0 writes list of global bodies to file
------------------------------------------------------------------------- */

void FixABF::write_restart_file(const char *file)
{
  if (comm->me) return;

  auto outfile = std::string(file) + ".rigid";
  FILE *fp = fopen(outfile.c_str(), "w");
  if (fp == nullptr)
    error->one(FLERR, "Cannot open fix rigid restart file {}: {}", outfile, utils::getsyserror());

  fmt::print(fp, "# fix rigid mass, COM, inertia tensor info for {} bodies on timestep {}\n\n",
             nbody, update->ntimestep);
  fmt::print(fp, "{}\n", nbody);

  // compute I tensor against xyz axes from diagonalized I and current quat
  // Ispace = P Idiag P_transpose
  // P is stored column-wise in exyz_space

  int xbox, ybox, zbox;
  double p[3][3], pdiag[3][3], ispace[3][3];

  int id;
  for (int i = 0; i < nbody; i++) {
    if (rstyle == SINGLE || rstyle == GROUP)
      id = i + 1;
    else
      id = body2mol[i];

    MathExtra::col2mat(ex_space[i], ey_space[i], ez_space[i], p);
    MathExtra::times3_diag(p, inertia[i], pdiag);
    MathExtra::times3_transpose(pdiag, p, ispace);

    xbox = (imagebody[i] & IMGMASK) - IMGMAX;
    ybox = (imagebody[i] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (imagebody[i] >> IMG2BITS) - IMGMAX;

    fprintf(fp,
            "%d %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e "
            "%-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %d %d %d\n",
            id, masstotal[i], xcm[i][0], xcm[i][1], xcm[i][2], ispace[0][0], ispace[1][1],
            ispace[2][2], ispace[0][1], ispace[0][2], ispace[1][2], vcm[i][0], vcm[i][1], vcm[i][2],
            angmom[i][0], angmom[i][1], angmom[i][2], xbox, ybox, zbox);
  }

  fclose(fp);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixABF::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = (double) nmax * sizeof(int);
  bytes += (double) nmax * sizeof(imageint);
  bytes += (double) nmax * 3 * sizeof(double);
  bytes += (double) maxvatom * 6 * sizeof(double);    // vatom
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixABF::grow_arrays(int nmax)
{
  memory->grow(body, nmax, "abf:body");
  memory->grow(xcmimage, nmax, "abf:xcmimage");
  memory->grow(displace, nmax, 3, "abf:displace");

  // check for regrow of vatom
  // must be done whether per-atom virial is accumulated on this step or not
  //   b/c this is only time grow_array() may be called
  // need to regrow b/c vatom is calculated before and after atom migration

  if (nmax > maxvatom) {
    maxvatom = atom->nmax;
    memory->grow(vatom, maxvatom, 6, "fix:vatom");
  }
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixABF::copy_arrays(int i, int j, int /*delflag*/)
{
  body[j] = body[i];
  xcmimage[j] = xcmimage[i];
  displace[j][0] = displace[i][0];
  displace[j][1] = displace[i][1];
  displace[j][2] = displace[i][2];

  // must also copy vatom if per-atom virial calculated on this timestep
  // since vatom is calculated before and after atom migration

  if (vflag_atom)
    for (int k = 0; k < 6; k++) vatom[j][k] = vatom[i][k];
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixABF::set_arrays(int i)
{
  body[i] = -1;
  xcmimage[i] = 0;
  displace[i][0] = 0.0;
  displace[i][1] = 0.0;
  displace[i][2] = 0.0;

  // must also zero vatom if per-atom virial calculated on this timestep
  // since vatom is calculated before and after atom migration

  if (vflag_atom)
    for (int k = 0; k < 6; k++) vatom[i][k] = 0.0;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixABF::pack_exchange(int i, double *buf)
{
  buf[0] = ubuf(body[i]).d;
  buf[1] = ubuf(xcmimage[i]).d;
  buf[2] = displace[i][0];
  buf[3] = displace[i][1];
  buf[4] = displace[i][2];

  // must also pack vatom if per-atom virial calculated on this timestep
  // since vatom is calculated before and after atom migration

  int m = 5;
  if (vflag_atom)
    for (int k = 0; k < 6; k++) { buf[m++] = vatom[i][k]; }

  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixABF::unpack_exchange(int nlocal, double *buf)
{
  body[nlocal] = (int) ubuf(buf[0]).i;
  xcmimage[nlocal] = (imageint) ubuf(buf[1]).i;
  displace[nlocal][0] = buf[2];
  displace[nlocal][1] = buf[3];
  displace[nlocal][2] = buf[4];

  // must also unpack vatom if per-atom virial calculated on this timestep
  // since vatom is calculated before and after atom migration

  int m = 5;
  if (vflag_atom)
    for (int k = 0; k < 6; k++) vatom[nlocal][k] = buf[m++];

  return m;
}

/* ---------------------------------------------------------------------- */

void FixABF::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  dtq = 0.5 * update->dt;
}

/* ----------------------------------------------------------------------
   zero linear momentum of each rigid body
   set Vcm to 0.0, then reset velocities of particles via set_v()
------------------------------------------------------------------------- */

void FixABF::zero_momentum()
{
  for (int ibody = 0; ibody < nbody; ibody++) vcm[ibody][0] = vcm[ibody][1] = vcm[ibody][2] = 0.0;

  evflag = 0;
  set_v();
}

/* ----------------------------------------------------------------------
   zero angular momentum of each rigid body
   set angmom/omega to 0.0, then reset velocities of particles via set_v()
------------------------------------------------------------------------- */

void FixABF::zero_rotation()
{
  for (int ibody = 0; ibody < nbody; ibody++) {
    angmom[ibody][0] = angmom[ibody][1] = angmom[ibody][2] = 0.0;
    omega[ibody][0] = omega[ibody][1] = omega[ibody][2] = 0.0;
  }

  evflag = 0;
  set_v();
}

/* ---------------------------------------------------------------------- */

int FixABF::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0], "bodyforces") == 0) {
    if (narg < 2) error->all(FLERR, "Illegal fix_modify command");
    if (strcmp(arg[1], "early") == 0)
      earlyflag = 1;
    else if (strcmp(arg[1], "late") == 0)
      earlyflag = 0;
    else
      error->all(FLERR, "Illegal fix_modify command");

    // reset fix mask
    // must do here and not in init,
    // since modify.cpp::init() uses fix masks before calling fix::init()

    for (int i = 0; i < modify->nfix; i++)
      if (strcmp(modify->fix[i]->id, id) == 0) {
        if (earlyflag) modify->fmask[i] |= POST_FORCE;
        break;
      }
    return 2;
  }

  return 0;
}

/* ----------------------------------------------------------------------
   return temperature of collection of rigid bodies
------------------------------------------------------------------------- */

double FixABF::compute_scalar()
{
  double wbody[3], rot[3][3];

  double t = 0.0;
  for (int i = 0; i < nbody; i++) {
    t += masstotal[i] * (vcm[i][0] * vcm[i][0] + vcm[i][1] * vcm[i][1] + vcm[i][2] * vcm[i][2]);

    // wbody = angular velocity in body frame

    MathExtra::quat_to_mat(quat[i], rot);
    MathExtra::transpose_matvec(rot, angmom[i], wbody);
    if (inertia[i][0] == 0.0)
      wbody[0] = 0.0;
    else
      wbody[0] /= inertia[i][0];
    if (inertia[i][1] == 0.0)
      wbody[1] = 0.0;
    else
      wbody[1] /= inertia[i][1];
    if (inertia[i][2] == 0.0)
      wbody[2] = 0.0;
    else
      wbody[2] /= inertia[i][2];

    t += inertia[i][0] * wbody[0] * wbody[0] + inertia[i][1] * wbody[1] * wbody[1] +
        inertia[i][2] * wbody[2] * wbody[2];
  }

  t *= tfactor;
  return t;
}

/* ---------------------------------------------------------------------- */

void *FixABF::extract(const char *str, int &dim)
{
  dim = 0;

  if (strcmp(str, "body") == 0) {
    if (!setupflag) return nullptr;
    dim = 1;
    return body;
  }
  if (strcmp(str, "masstotal") == 0) {
    if (!setupflag) return nullptr;
    dim = 1;
    return masstotal;
  }
  if (strcmp(str, "t_target") == 0) {
    dim = 0;
    return &t_target;
  }

  return nullptr;
}

/* ----------------------------------------------------------------------
   return translational KE for all rigid bodies
   KE = 1/2 M Vcm^2
------------------------------------------------------------------------- */

double FixABF::extract_ke()
{
  double ke = 0.0;
  for (int i = 0; i < nbody; i++)
    ke += masstotal[i] * (vcm[i][0] * vcm[i][0] + vcm[i][1] * vcm[i][1] + vcm[i][2] * vcm[i][2]);

  return 0.5 * ke;
}

/* ----------------------------------------------------------------------
   return rotational KE for all rigid bodies
   Erotational = 1/2 I wbody^2
------------------------------------------------------------------------- */

double FixABF::extract_erotational()
{
  double wbody[3], rot[3][3];

  double erotate = 0.0;
  for (int i = 0; i < nbody; i++) {

    // wbody = angular velocity in body frame

    MathExtra::quat_to_mat(quat[i], rot);
    MathExtra::transpose_matvec(rot, angmom[i], wbody);
    if (inertia[i][0] == 0.0)
      wbody[0] = 0.0;
    else
      wbody[0] /= inertia[i][0];
    if (inertia[i][1] == 0.0)
      wbody[1] = 0.0;
    else
      wbody[1] /= inertia[i][1];
    if (inertia[i][2] == 0.0)
      wbody[2] = 0.0;
    else
      wbody[2] /= inertia[i][2];

    erotate += inertia[i][0] * wbody[0] * wbody[0] + inertia[i][1] * wbody[1] * wbody[1] +
        inertia[i][2] * wbody[2] * wbody[2];
  }

  return 0.5 * erotate;
}

/* ----------------------------------------------------------------------
   return attributes of a rigid body
   15 values per body
   xcm = 0,1,2; vcm = 3,4,5; fcm = 6,7,8; torque = 9,10,11; image = 12,13,14
------------------------------------------------------------------------- */

double FixABF::compute_array(int i, int j)
{
  if (j < 3) return xcm[i][j];
  if (j < 6) return vcm[i][j - 3];
  if (j < 9) return fcm[i][j - 6];
  if (j < 12) return torque[i][j - 9];
  if (j == 12) return (imagebody[i] & IMGMASK) - IMGMAX;
  if (j == 13) return (imagebody[i] >> IMGBITS & IMGMASK) - IMGMAX;
  return (imagebody[i] >> IMG2BITS) - IMGMAX;
}

void FixABF::broadcastState()
{
  int rank;
  MPI_Comm_rank(abfComm, &rank);

  int root = 0;
  if (isRankCloseToAbf && rank == 0) root = me;

  MPI_Allreduce(MPI_IN_PLACE, &root, 1, MPI_INT, MPI_SUM, world);

  MPI_Bcast(xcm[0], 3, MPI_DOUBLE, root, world);
  MPI_Bcast(vcm[0], 3, MPI_DOUBLE, root, world);
  MPI_Bcast(fcm[0], 3, MPI_DOUBLE, root, world);
  MPI_Bcast(ex_space[0], 3, MPI_DOUBLE, root, world);
  MPI_Bcast(ey_space[0], 3, MPI_DOUBLE, root, world);
  MPI_Bcast(ez_space[0], 3, MPI_DOUBLE, root, world);
  MPI_Bcast(angmom[0], 3, MPI_DOUBLE, root, world);
  MPI_Bcast(omega[0], 3, MPI_DOUBLE, root, world);
  MPI_Bcast(quat[0], 4, MPI_DOUBLE, root, world);
}

static bool AABB(const double *loa, const double *hia, const double *lob, const double *hib)
{
  return loa[0] <= hib[0] && hia[0] >= lob[0] && loa[1] <= hib[1] && hia[1] >= lob[1] &&
      loa[2] <= hib[2] && hia[2] >= lob[2];
}

void FixABF::updateAbfComm()
{

  if (nbody != 1) error->all(FLERR, "fix abf requires exactly one object");
  // get extents of the ABF

  const int nlocal = atom->nlocal;
  double **x = atom->x;

  double drmin[3] = {1e9, 1e9, 1e9};
  double drmax[3] = {-1e9, -1e9, -1e9};

  for (int i = 0; i < nlocal; i++) {
    if (body[i] < 0) continue;
    const int ibody = body[i];

    double unwrap[3];
    domain->unmap(x[i], xcmimage[i], unwrap);

    for (int j = 0; j < 3; ++j) {
      const double dj = unwrap[j] - xcm[ibody][j];

      drmin[j] = std::min(drmin[j], dj);
      drmax[j] = std::max(drmax[j], dj);
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, drmin, 3, MPI_DOUBLE, MPI_MIN, world);
  MPI_Allreduce(MPI_IN_PLACE, drmax, 3, MPI_DOUBLE, MPI_MAX, world);

  double cm[3] = {xcm[0][0], xcm[0][1], xcm[0][2]};
  // domain->minimum_image(cm[0], cm[1], cm[2]);

  double loAbf[3], hiAbf[3];
  for (int j = 0; j < 3; ++j) {
    loAbf[j] = cm[j] + drmin[j] - abfCommMargin;
    hiAbf[j] = cm[j] + drmax[j] + abfCommMargin;
  }

  isRankCloseToAbf = AABB(loAbf, hiAbf, domain->sublo, domain->subhi);

  const int color = isRankCloseToAbf ? 0 : 1;
  const int key = me;

  // printf("%d %d   %g %g %g   %g %g %g -  %g %g %g   %g %g %g\n", me, color,
  //        loAbf[0], loAbf[1], loAbf[2], hiAbf[0], hiAbf[1], hiAbf[2],
  //        domain->sublo[0], domain->sublo[1], domain->sublo[2],
  //        domain->subhi[0], domain->subhi[1], domain->subhi[2]);

  if (abfComm != MPI_COMM_NULL) MPI_Comm_free(&abfComm);

  MPI_Comm_split(world, color, key, &abfComm);
}
