// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "set.h"

#include "arg_info.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "random_mars.h"
#include "random_park.h"
#include "region.h"
#include "variable.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

enum{ATOM_SELECT,TYPE_SELECT,GROUP_SELECT,REGION_SELECT};

enum{TYPE,TYPE_FRACTION,TYPE_RATIO,TYPE_SUBSET,
     MOLECULE,X,Y,Z,VX,VY,VZ,CHARGE,MASS,SHAPE,LENGTH,TRI,
     DIPOLE,DIPOLE_RANDOM,SPIN_ATOM,SPIN_RANDOM,SPIN_ELECTRON,RADIUS_ELECTRON,
     QUAT,QUAT_RANDOM,THETA,THETA_RANDOM,ANGMOM,OMEGA,TEMPERATURE,
     DIAMETER,RADIUS_ATOM,DENSITY,VOLUME,IMAGE,BOND,
     SPH_E,SPH_CV,SPH_RHO,EDPD_TEMP,EDPD_CV,CC,SMD_MASS_DENSITY,
     SMD_CONTACT_RADIUS,DPDTHETA,EPSILON,IVEC,DVEC,IARRAY,DARRAY};

#define BIG INT_MAX

/* ---------------------------------------------------------------------- */

void Set::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Set command before simulation box is defined");
  if (atom->natoms == 0)
    error->all(FLERR,"Set command on system without atoms");
  if (narg < 4) error->all(FLERR,"Illegal set command: need at least four arguments");

  // style and ID info

  if (strcmp(arg[0],"atom") == 0) style = ATOM_SELECT;
  else if (strcmp(arg[0],"type") == 0) style = TYPE_SELECT;
  else if (strcmp(arg[0],"group") == 0) style = GROUP_SELECT;
  else if (strcmp(arg[0],"region") == 0) style = REGION_SELECT;
  else error->all(FLERR,"Unknown set command style: {}", arg[0]);

  id = utils::strdup(arg[1]);
  select = nullptr;
  selection(atom->nlocal);

  // loop over keyword/value pairs
  // call appropriate routine to reset attributes

  if (comm->me == 0) utils::logmesg(lmp,"Setting atom values ...\n");

  int allcount,origarg;

  int iarg = 2;
  while (iarg < narg) {
    varflag = varflag1 = varflag2 = varflag3 = varflag4 = 0;
    count = 0;
    origarg = iarg;

    if (strcmp(arg[iarg],"mass") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set mass", error);
      if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
      else dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (!atom->rmass_flag)
        error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
      set(MASS);
      iarg += 2;
    } else if (strcmp(arg[iarg],"density") == 0 ||
               (strcmp(arg[iarg],"density/disc") == 0)) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set density", error);
      if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
      else dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (!atom->rmass_flag)
        error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
      if (dvalue <= 0.0) error->all(FLERR,"Invalid density in set command");
      discflag = 0;
      if (strcmp(arg[iarg],"density/disc") == 0) {
        discflag = 1;
        if (domain->dimension != 2)
          error->all(FLERR,"Density/disc option requires 2d simulation");
      }
      set(DENSITY);
      iarg += 2;

    } else if (strcmp(arg[iarg],"image") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR, "set image", error);
      ximageflag = yimageflag = zimageflag = 0;
      if (strcmp(arg[iarg+1],"NULL") != 0) {
        ximageflag = 1;
        if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
        else ximage = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      }
      if (strcmp(arg[iarg+2],"NULL") != 0) {
        yimageflag = 1;
        if (utils::strmatch(arg[iarg+2],"^v_")) varparse(arg[iarg+2],2);
        else yimage = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
      }
      if (strcmp(arg[iarg+3],"NULL") != 0) {
        zimageflag = 1;
        if (utils::strmatch(arg[iarg+3],"^v_")) varparse(arg[iarg+3],3);
        else zimage = utils::inumeric(FLERR,arg[iarg+3],false,lmp);
      }
      if (ximageflag && ximage && !domain->xperiodic)
        error->all(FLERR,
                   "Cannot set non-zero image flag for non-periodic dimension");
      if (yimageflag && yimage && !domain->yperiodic)
        error->all(FLERR,
                   "Cannot set non-zero image flag for non-periodic dimension");
      if (zimageflag && zimage && !domain->zperiodic)
        error->all(FLERR,
                   "Cannot set non-zero image flag for non-periodic dimension");
      set(IMAGE);
      iarg += 4;

    } else if (strcmp(arg[iarg],"sph/e") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set sph/e", error);
      if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
      else dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (!atom->esph_flag)
        error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
      set(SPH_E);
      iarg += 2;

    } else if (strcmp(arg[iarg],"sph/cv") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set sph/cv", error);
      if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
      else dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (!atom->cv_flag)
        error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
      set(SPH_CV);
      iarg += 2;

    } else if (strcmp(arg[iarg],"sph/rho") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set sph/rho", error);
      if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
      else dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (!atom->rho_flag)
        error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
      set(SPH_RHO);
      iarg += 2;

    } else if (strcmp(arg[iarg],"edpd/temp") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set edpd/temp", error);
      if (strcmp(arg[iarg+1],"NULL") == 0) dvalue = -1.0;
      else if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
      else {
        dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        if (dvalue < 0.0) error->all(FLERR,"Illegal set command");
      }
      if (!atom->edpd_flag)
        error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
      set(EDPD_TEMP);
      iarg += 2;

    } else if (strcmp(arg[iarg],"edpd/cv") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set edpd/cv", error);
      if (strcmp(arg[iarg+1],"NULL") == 0) dvalue = -1.0;
      else if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
      else {
        dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        if (dvalue < 0.0) error->all(FLERR,"Illegal set command");
      }
      if (!atom->edpd_flag)
        error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
      set(EDPD_CV);
      iarg += 2;

    } else if (strcmp(arg[iarg],"cc") == 0) {
      if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "set cc", error);
      if (strcmp(arg[iarg+1],"NULL") == 0) dvalue = -1.0;
      else if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
      else {
        cc_index = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
        dvalue = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        if (cc_index < 1) error->all(FLERR,"Illegal set command");
      }
      if (!atom->tdpd_flag)
        error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
      set(CC);
      iarg += 3;

    } else if (strcmp(arg[iarg],"dpd/theta") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set dpd/theta", error);
      if (strcmp(arg[iarg+1],"NULL") == 0) dvalue = -1.0;
      else if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
      else {
        dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        if (dvalue < 0.0) error->all(FLERR,"Illegal set command");
      }
      if (!atom->dpd_flag)
        error->all(FLERR,"Cannot set attribute {} for atom style {}", arg[iarg], atom->get_style());
      set(DPDTHETA);
      iarg += 2;

    } else {

      // set custom per-atom vector or array or error out

      int flag,cols;
      ArgInfo argi(arg[iarg],ArgInfo::DNAME|ArgInfo::INAME);
      const char *pname = argi.get_name();
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "set", error);
      index_custom = atom->find_custom(argi.get_name(),flag,cols);
      if (index_custom < 0)
        error->all(FLERR,"Set keyword or custom property {} does not exist",pname);

      switch (argi.get_type()) {

      case ArgInfo::INAME:
        if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
        else ivalue = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
        if (flag != 0) error->all(FLERR,"Set command custom property {} is not integer",pname);

        if (argi.get_dim() == 0) {
          if (cols > 0)
            error->all(FLERR,"Set command custom integer property {} is not a vector",pname);
          set(IVEC);
        } else if (argi.get_dim() == 1) {
          if (cols == 0)
            error->all(FLERR,"Set command custom integer property {} is not an array",pname);
          icol_custom = argi.get_index1();
          if (icol_custom <= 0 || icol_custom > cols)
            error->all(FLERR,"Set command per-atom custom integer array {} is accessed "
                       "out-of-range",pname);
          set(IARRAY);
        } else error->all(FLERR,"Illegal set command");
        break;

      case ArgInfo::DNAME:
        if (utils::strmatch(arg[iarg+1],"^v_")) varparse(arg[iarg+1],1);
        else dvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        if (flag != 1) error->all(FLERR,"Custom property {} is not floating-point",argi.get_name());

        if (argi.get_dim() == 0) {
          if (cols > 0)
            error->all(FLERR,"Set command custom double property {} is not a vector",pname);
          set(DVEC);
        } else if (argi.get_dim() == 1) {
          if (cols == 0)
            error->all(FLERR,"Set command custom double property {} is not an array",pname);
          icol_custom = argi.get_index1();
          if (icol_custom <= 0 || icol_custom > cols)
            error->all(FLERR,"Set command per-atom custom double array {} is "
                       "accessed out-of-range",pname);
          set(DARRAY);
        } else error->all(FLERR,"Illegal set command");
        break;

      default:
        error->all(FLERR,"Illegal set command");
        break;
      }
      iarg += 2;
    }

    // statistics
    // for CC option, include species index

    MPI_Allreduce(&count,&allcount,1,MPI_INT,MPI_SUM,world);

    if (comm->me == 0) {
      if (strcmp(arg[origarg],"cc") == 0)
        utils::logmesg(lmp,"  {} settings made for {} index {}\n",
                       allcount,arg[origarg],arg[origarg+1]);
      else
        utils::logmesg(lmp,"  {} settings made for {}\n",
                       allcount,arg[origarg]);
    }
  }

  // free local memory

  delete[] id;
  delete[] select;
}

/* ----------------------------------------------------------------------
   select atoms according to ATOM, MOLECULE, TYPE, GROUP, REGION style
   n = nlocal or nlocal+nghost depending on keyword
------------------------------------------------------------------------- */

void Set::selection(int n)
{
  delete[] select;
  select = new int[n];
  int nlo,nhi;

  if (style == ATOM_SELECT) {
    if (atom->tag_enable == 0)
      error->all(FLERR,"Cannot use set atom with no atom IDs defined");
    bigint nlobig,nhibig;
    utils::bounds(FLERR,id,1,MAXTAGINT,nlobig,nhibig,error);

    tagint *tag = atom->tag;
    for (int i = 0; i < n; i++)
      if (tag[i] >= nlobig && tag[i] <= nhibig) select[i] = 1;
      else select[i] = 0;

  } else if (style == TYPE_SELECT) {
    int *type = atom->type;
    for (int i = 0; i < n; i++)
      if (type[i] >= nlo && type[i] <= nhi) select[i] = 1;
      else select[i] = 0;

  } else if (style == GROUP_SELECT) {
    int igroup = group->find(id);
    if (igroup == -1) error->all(FLERR,"Could not find set group ID {}", id);
    int groupbit = group->bitmask[igroup];

    int *mask = atom->mask;
    for (int i = 0; i < n; i++)
      if (mask[i] & groupbit) select[i] = 1;
      else select[i] = 0;

  } else if (style == REGION_SELECT) {
    auto region = domain->get_region_by_id(id);
    if (!region) error->all(FLERR,"Set region {} does not exist", id);
    region->prematch();

    double **x = atom->x;
    for (int i = 0; i < n; i++)
      if (region->match(x[i][0],x[i][1],x[i][2]))
        select[i] = 1;
      else select[i] = 0;
  }
}

/* ----------------------------------------------------------------------
   set owned atom properties directly
   either scalar or per-atom values from atom-style variable(s)
------------------------------------------------------------------------- */

void Set::set(int keyword)
{
  // evaluate atom-style variable(s) if necessary

  vec1 = vec2 = vec3 = vec4 = nullptr;

  if (varflag) {
    int nlocal = atom->nlocal;
    if (varflag1) {
      memory->create(vec1,nlocal,"set:vec1");
      input->variable->compute_atom(ivar1,0,vec1,1,0);
    }
    if (varflag2) {
      memory->create(vec2,nlocal,"set:vec2");
      input->variable->compute_atom(ivar2,0,vec2,1,0);
    }
    if (varflag3) {
      memory->create(vec3,nlocal,"set:vec3");
      input->variable->compute_atom(ivar3,0,vec3,1,0);
    }
    if (varflag4) {
      memory->create(vec4,nlocal,"set:vec4");
      input->variable->compute_atom(ivar4,0,vec4,1,0);
    }
  }

  // check if properties of atoms in rigid bodies are updated
  // that are cached as per-body data.
  switch (keyword) {
  case X:
  case Y:
  case Z:
  case MOLECULE:
  case MASS:
  case ANGMOM:
  case SHAPE:
  case DIAMETER:
  case DENSITY:
  case TEMPERATURE:
  case QUAT:
  case IMAGE:
    if (modify->check_rigid_list_overlap(select))
      error->warning(FLERR,"Changing a property of atoms in rigid bodies "
                     "that has no effect unless rigid bodies are rebuild");
    break;
  default: // assume no conflict for all other properties
    break;
  }

  // loop over selected atoms

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    if (!select[i]) continue;

    // overwrite dvalue, ivalue, xyzw value if variables defined
    // else the input script scalar value remains in place

    if (varflag) {
      if (varflag1) {
        dvalue = xvalue = vec1[i];
        ivalue = static_cast<int> (dvalue);
      }
      if (varflag2) yvalue = vec2[i];
      if (varflag3) zvalue = vec3[i];
      if (varflag4) wvalue = vec4[i];
    }

    // set values in per-atom arrays
    // error check here in case atom-style variables generated bogus value

    if (keyword == TYPE) {
      if (ivalue <= 0 || ivalue > atom->ntypes)
        error->one(FLERR,"Invalid value in set command");
      atom->type[i] = ivalue;
    }
    else if (keyword == X) atom->x[i][0] = dvalue;
    else if (keyword == Y) atom->x[i][1] = dvalue;
    else if (keyword == Z) atom->x[i][2] = dvalue;
    else if (keyword == VX) atom->v[i][0] = dvalue;
    else if (keyword == VY) atom->v[i][1] = dvalue;
    else if (keyword == VZ) atom->v[i][2] = dvalue;
    else if (keyword == EDPD_TEMP) atom->edpd_temp[i] = dvalue;
    else if (keyword == EDPD_CV) atom->edpd_cv[i] = dvalue;
    else if (keyword == CC) atom->cc[i][cc_index-1] = dvalue;

    else if (keyword == SMD_MASS_DENSITY) {
      // set mass from volume and supplied mass density
      atom->rmass[i] = atom->vfrac[i] * dvalue;
    }
    else if (keyword == SMD_CONTACT_RADIUS) atom->contact_radius[i] = dvalue;

    else if (keyword == IMAGE) {
      int xbox = (atom->image[i] & IMGMASK) - IMGMAX;
      int ybox = (atom->image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      int zbox = (atom->image[i] >> IMG2BITS) - IMGMAX;
      if (varflag1) ximage = static_cast<int>(xvalue);
      if (varflag2) yimage = static_cast<int>(yvalue);
      if (varflag3) zimage = static_cast<int>(zvalue);
      if (ximageflag) xbox = ximage;
      if (yimageflag) ybox = yimage;
      if (zimageflag) zbox = zimage;
      atom->image[i] = ((imageint) (xbox + IMGMAX) & IMGMASK) |
        (((imageint) (ybox + IMGMAX) & IMGMASK) << IMGBITS) |
        (((imageint) (zbox + IMGMAX) & IMGMASK) << IMG2BITS);
    }
    count++;
  }
  // clear up per-atom memory if allocated

  memory->destroy(vec1);
  memory->destroy(vec2);
  memory->destroy(vec3);
  memory->destroy(vec4);
}

/* ---------------------------------------------------------------------- */

void Set::varparse(const char *name, int m)
{
  varflag = 1;
  int ivar = input->variable->find(name+2);

  if (ivar < 0)
    error->all(FLERR,"Variable name {} for set command does not exist", name);
  if (!input->variable->atomstyle(ivar))
    error->all(FLERR,"Variable {} for set command is invalid style", name);

  if (m == 1) {
    varflag1 = 1; ivar1 = ivar;
  } else if (m == 2) {
    varflag2 = 1; ivar2 = ivar;
  } else if (m == 3) {
    varflag3 = 1; ivar3 = ivar;
  } else if (m == 4) {
    varflag4 = 1; ivar4 = ivar;
  }
}
