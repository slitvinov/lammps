#include "run.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "integrate.h"
#include "modify.h"
#include "update.h"
#include <cstring>
using namespace LAMMPS_NS;
Run::Run(LAMMPS *lmp) : Command(lmp) {}
void Run::command(int narg, char **arg)
{
  if (narg < 1) utils::missing_cmd_args(FLERR, "run", error);
  if (domain->box_exist == 0)
    error->all(FLERR,"Run command before simulation box is defined");
  bigint nsteps_input = utils::bnumeric(FLERR,arg[0],false,lmp);
  int uptoflag = 0;
  int startflag = 0;
  int stopflag = 0;
  bigint start,stop;
  int preflag = 1;
  int postflag = 1;
  int nevery = 0;
  int ncommands = 0;
  int first,last;
  int iarg = 1;
  int nsteps;
  if (!uptoflag) {
    if (nsteps_input < 0 || nsteps_input > MAXSMALLINT)
      error->all(FLERR,"Invalid run command N value: {}", nsteps_input);
    nsteps = static_cast<int> (nsteps_input);
  }
  char **commands = nullptr;
  update->whichflag = 1;
  if (nevery == 0) {
    update->nsteps = nsteps;
    update->firststep = update->ntimestep;
    update->laststep = update->ntimestep + nsteps;
    if (update->laststep < 0 || update->laststep < update->firststep)
      error->all(FLERR,"Too many timesteps");
    if (startflag) update->beginstep = start;
    else update->beginstep = update->firststep;
    if (stopflag) update->endstep = stop;
    else update->endstep = update->laststep;
    if (preflag || update->first_update == 0) {
      lmp->init();
      update->integrate->setup(1);
    }
    update->integrate->run(nsteps);
    update->integrate->cleanup();
  }
  update->whichflag = 0;
  update->firststep = update->laststep = 0;
  update->beginstep = update->endstep = 0;
  if (commands) {
    for (int i = 0; i < ncommands; i++) delete [] commands[i];
    delete [] commands;
  }
}
