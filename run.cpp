#include <cstring>
#include <unordered_set>
#include <map>
#include <vector>
#include <string>
#include <cmath>
#include "utils.h"
#include "lammps.h"
#include "pointers.h"
#include "command.h"
#include "run.h"
#include "domain.h"
#include "input.h"
#include "integrate.h"
#include "modify.h"
#include "update.h"
using namespace LAMMPS_NS;
Run::Run(LAMMPS *lmp) : Command(lmp) {}
void Run::command(int narg, char **arg) {
  bigint nsteps_input = utils::bnumeric(FLERR, arg[0], false, lmp);
  bigint start, stop;
  int preflag = 1;
  int postflag = 1;
  int nevery = 0;
  int ncommands = 0;
  int first, last;
  int iarg = 1;
  int nsteps;
  nsteps = static_cast<int>(nsteps_input);
  update->whichflag = 1;
  update->nsteps = nsteps;
  update->firststep = update->ntimestep;
  update->laststep = update->ntimestep + nsteps;
  update->beginstep = update->firststep;
  update->endstep = update->laststep;
  if (preflag || update->first_update == 0) {
    lmp->init();
    update->integrate->setup(1);
  }
  update->integrate->run(nsteps);
  update->integrate->cleanup();
  update->whichflag = 0;
  update->firststep = update->laststep = 0;
  update->beginstep = update->endstep = 0;
}
