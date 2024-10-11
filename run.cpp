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
  while (iarg < narg) {
    if (strcmp(arg[iarg],"upto") == 0) {
      if (iarg+1 > narg) utils::missing_cmd_args(FLERR, "run upto", error);
      uptoflag = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg],"start") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "run start", error);
      startflag = 1;
      start = utils::bnumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"stop") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "run stop", error);
      stopflag = 1;
      stop = utils::bnumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"pre") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "run pre", error);
      preflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"post") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "run post", error);
      postflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "run every", error);
      nevery = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (nevery <= 0) error->all(FLERR, "Invalid run every argument: {}", nevery);
      first = iarg+2;
      last = narg-1;
      ncommands = last-first + 1;
      if (ncommands == 1 && strcmp(arg[first],"NULL") == 0) ncommands = 0;
      iarg = narg;
    } else error->all(FLERR,"Unknown run keyword: {}", arg[iarg]);
  }
  int nsteps;
  if (!uptoflag) {
    if (nsteps_input < 0 || nsteps_input > MAXSMALLINT)
      error->all(FLERR,"Invalid run command N value: {}", nsteps_input);
    nsteps = static_cast<int> (nsteps_input);
  } else {
    bigint delta = nsteps_input - update->ntimestep;
    if (delta < 0 || delta > MAXSMALLINT)
      error->all(FLERR,"Invalid run command upto value: {}", delta);
    nsteps = static_cast<int> (delta);
  }
  if (startflag) {
    if (start < 0)
      error->all(FLERR,"Invalid run command start value: {}", start);
    if (start > update->ntimestep)
      error->all(FLERR,"Run command start value is after start of run");
  }
  if (stopflag) {
    if (stop < 0)
      error->all(FLERR,"Invalid run command stop value: {}", stop);
    if (stop < update->ntimestep + nsteps)
      error->all(FLERR,"Run command stop value is before end of run");
  }
  if (!preflag && utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Run flag 'pre no' not compatible with r-RESPA");
  char **commands = nullptr;
  if (nevery && ncommands > 0) {
    commands = new char*[ncommands];
    ncommands = 0;
    for (int i = first; i <= last; i++) {
      commands[ncommands] = utils::strdup(arg[i]);
      ncommands++;
    }
  }
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
  } else {
    int iter = 0;
    int nleft = nsteps;
    while (nleft > 0 || iter == 0) {
      nsteps = MIN(nleft,nevery);
      update->nsteps = nsteps;
      update->firststep = update->ntimestep;
      update->laststep = update->ntimestep + nsteps;
      if (update->laststep < 0 || update->laststep < update->firststep)
        error->all(FLERR,"Too many timesteps");
      if (startflag) update->beginstep = start;
      else update->beginstep = update->firststep;
      if (stopflag) update->endstep = stop;
      else update->endstep = update->laststep;
      if (preflag || iter == 0) {
        lmp->init();
        update->integrate->setup(1);
      }
      update->integrate->run(nsteps);
      update->integrate->cleanup();
      if (ncommands) {
        modify->clearstep_compute();
        for (int i = 0; i < ncommands; i++) input->one(commands[i]);
        modify->addstep_compute(update->ntimestep + nevery);
      }
      nleft -= nsteps;
      iter++;
    }
  }
  update->whichflag = 0;
  update->firststep = update->laststep = 0;
  update->beginstep = update->endstep = 0;
  if (commands) {
    for (int i = 0; i < ncommands; i++) delete [] commands[i];
    delete [] commands;
  }
}
