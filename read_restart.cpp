#include "read_restart.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_read_restart.h"
#include "force.h"
#include "group.h"
#include "label_map.h"
#include "memory.h"
#include "modify.h"
#include "mpiio.h"
#include "pair.h"
#include "update.h"
#include <cstring>
#include "lmprestart.h"
using namespace LAMMPS_NS;
ReadRestart::ReadRestart(LAMMPS *lmp) : Command(lmp), mpiio(nullptr) {}
void ReadRestart::command(int narg, char **arg)
{
  if (narg != 1 && narg != 2) error->all(FLERR,"Illegal read_restart command");
  if (domain->box_exist)
    error->all(FLERR,"Cannot read_restart after simulation box is defined");
  MPI_Barrier(world);
  double time1 = platform::walltime();
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  int remapflag = 1;
  if (narg == 2) {
    if (strcmp(arg[1],"noremap") == 0) remapflag = 0;
    else if (strcmp(arg[1],"remap") == 0) remapflag = 1;
    else error->all(FLERR,"Illegal read_restart command");
  }
  char *file;
  if (strchr(arg[0],'*')) {
    int n=0;
    if (me == 0) {
      auto fn = file_search(arg[0]);
      n = fn.size()+1;
      file = utils::strdup(fn);
    }
    MPI_Bcast(&n,1,MPI_INT,0,world);
    if (me != 0) file = new char[n];
    MPI_Bcast(file,n,MPI_CHAR,0,world);
  } else file = utils::strdup(arg[0]);
  if (strchr(arg[0],'%')) multiproc = 1;
  else multiproc = 0;
  if (strstr(arg[0],".mpiio")) mpiioflag = 1;
  else mpiioflag = 0;
  if (multiproc && mpiioflag)
    error->all(FLERR,"Read restart MPI-IO input not allowed with % in filename");
  if (mpiioflag) {
    mpiio = new RestartMPIIO(lmp);
    if (!mpiio->mpiio_exists)
      error->all(FLERR,"Reading from MPI-IO filename when MPIIO package is not installed");
  }
  if (me == 0) {
    utils::logmesg(lmp,"Reading restart file ...\n");
    std::string hfile = file;
    if (multiproc) {
      hfile.replace(hfile.find('%'),1,"base");
    }
    fp = fopen(hfile.c_str(),"rb");
    if (fp == nullptr)
      error->one(FLERR,"Cannot open restart file {}: {}", hfile, utils::getsyserror());
  }
  magic_string();
  endian();
  format_revision();
  check_eof_magic();
  if ((comm->me == 0) && (modify->get_fix_by_style("property/atom").size() > 0))
    error->warning(FLERR, "Fix property/atom command must be specified after read_restart "
                   "to restore its data.");
  header();
  domain->box_exist = 1;
  int n;
  if (nprocs == 1) n = static_cast<int> (atom->natoms);
  else n = static_cast<int> (LB_FACTOR * atom->natoms / nprocs);
  atom->allocate_type_arrays();
  bigint nbig = n;
  nbig = atom->avec->roundup(nbig);
  n = static_cast<int> (nbig);
  atom->avec->grow(n);
  n = atom->nmax;
  domain->print_box("  ");
  domain->set_initial_box(0);
  domain->set_global_box();
  comm->set_proc_grid();
  domain->set_local_box();
  group->read_restart(fp);
  type_arrays();
  int nextra = modify->read_restart(fp);
  atom->nextra_store = nextra;
  memory->create(atom->extra,n,nextra,"atom:extra");
  file_layout();
  if (multiproc && me == 0) {
    fclose(fp);
    fp = nullptr;
  }
  AtomVec *avec = atom->avec;
  int maxbuf = 0;
  double *buf = nullptr;
  int m,flag;
  if (mpiioflag) {
    mpiio->openForRead(file);
    memory->create(buf,assignedChunkSize,"read_restart:buf");
    mpiio->read((headerOffset+assignedChunkOffset),assignedChunkSize,buf);
    mpiio->close();
    if (!nextra) {
      atom->nlocal = 1;
      int perAtomSize = avec->size_restart();
      atom->nlocal = 0;
      int atomCt = (int) (assignedChunkSize / perAtomSize);
      if (atomCt > atom->nmax) avec->grow(atomCt);
    }
    m = 0;
    while (m < assignedChunkSize) m += avec->unpack_restart(&buf[m]);
  }
  else if (multiproc == 0) {
    int triclinic = domain->triclinic;
    imageint *iptr;
    double *x,lamda[3];
    double *coord,*sublo,*subhi;
    if (triclinic == 0) {
      sublo = domain->sublo;
      subhi = domain->subhi;
    } else {
      sublo = domain->sublo_lamda;
      subhi = domain->subhi_lamda;
    }
    for (int iproc = 0; iproc < nprocs_file; iproc++) {
      if (read_int() != PERPROC)
        error->all(FLERR,"Invalid flag in peratom section of restart file");
      n = read_int();
      if (n > maxbuf) {
        maxbuf = n;
        memory->destroy(buf);
        memory->create(buf,maxbuf,"read_restart:buf");
      }
      read_double_vec(n,buf);
      m = 0;
      while (m < n) {
        x = &buf[m+1];
        if (remapflag) {
          iptr = (imageint *) &buf[m+7];
          domain->remap(x,*iptr);
        }
        if (triclinic) {
          domain->x2lamda(x,lamda);
          coord = lamda;
        } else coord = x;
        if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
            coord[1] >= sublo[1] && coord[1] < subhi[1] &&
            coord[2] >= sublo[2] && coord[2] < subhi[2]) {
          m += avec->unpack_restart(&buf[m]);
        } else m += static_cast<int> (buf[m]);
      }
    }
    if (me == 0) {
      fclose(fp);
      fp = nullptr;
    }
  }
  else if (nprocs <= multiproc_file) {
    for (int iproc = me; iproc < multiproc_file; iproc += nprocs) {
      std::string procfile = file;
      procfile.replace(procfile.find('%'),1,fmt::format("{}",iproc));
      fp = fopen(procfile.c_str(),"rb");
      if (fp == nullptr)
        error->one(FLERR,"Cannot open restart file {}: {}",
                                     procfile, utils::getsyserror());
      utils::sfread(FLERR,&flag,sizeof(int),1,fp,nullptr,error);
      if (flag != PROCSPERFILE)
        error->one(FLERR,"Invalid flag in peratom section of restart file");
      int procsperfile;
      utils::sfread(FLERR,&procsperfile,sizeof(int),1,fp,nullptr,error);
      for (int i = 0; i < procsperfile; i++) {
        utils::sfread(FLERR,&flag,sizeof(int),1,fp,nullptr,error);
        if (flag != PERPROC)
          error->one(FLERR,"Invalid flag in peratom section of restart file");
        utils::sfread(FLERR,&n,sizeof(int),1,fp,nullptr,error);
        if (n > maxbuf) {
          maxbuf = n;
          memory->destroy(buf);
          memory->create(buf,maxbuf,"read_restart:buf");
        }
        utils::sfread(FLERR,buf,sizeof(double),n,fp,nullptr,error);
        m = 0;
        while (m < n) m += avec->unpack_restart(&buf[m]);
      }
      fclose(fp);
      fp = nullptr;
    }
  }
  else {
    int nfile = multiproc_file;
    int icluster = static_cast<int> ((bigint) me * nfile/nprocs);
    int fileproc = static_cast<int> ((bigint) icluster * nprocs/nfile);
    int fcluster = static_cast<int> ((bigint) fileproc * nfile/nprocs);
    if (fcluster < icluster) fileproc++;
    int fileprocnext =
      static_cast<int> ((bigint) (icluster+1) * nprocs/nfile);
    fcluster = static_cast<int> ((bigint) fileprocnext * nfile/nprocs);
    if (fcluster < icluster+1) fileprocnext++;
    int nclusterprocs = fileprocnext - fileproc;
    int filereader = 0;
    if (me == fileproc) filereader = 1;
    MPI_Comm clustercomm;
    MPI_Comm_split(world,icluster,0,&clustercomm);
    if (filereader) {
      std::string procfile = file;
      procfile.replace(procfile.find('%'),1,fmt::format("{}",icluster));
      fp = fopen(procfile.c_str(),"rb");
      if (fp == nullptr)
        error->one(FLERR,"Cannot open restart file {}: {}", procfile, utils::getsyserror());
    }
    int procsperfile;
    if (filereader) {
      utils::sfread(FLERR,&flag,sizeof(int),1,fp,nullptr,error);
      if (flag != PROCSPERFILE)
        error->one(FLERR,"Invalid flag in peratom section of restart file");
      utils::sfread(FLERR,&procsperfile,sizeof(int),1,fp,nullptr,error);
    }
    MPI_Bcast(&procsperfile,1,MPI_INT,0,clustercomm);
    int tmp,iproc;
    MPI_Request request;
    for (int i = 0; i < procsperfile; i++) {
      if (filereader) {
        utils::sfread(FLERR,&flag,sizeof(int),1,fp,nullptr,error);
        if (flag != PERPROC)
          error->one(FLERR,"Invalid flag in peratom section of restart file");
        utils::sfread(FLERR,&n,sizeof(int),1,fp,nullptr,error);
        if (n > maxbuf) {
          maxbuf = n;
          memory->destroy(buf);
          memory->create(buf,maxbuf,"read_restart:buf");
        }
        utils::sfread(FLERR,buf,sizeof(double),n,fp,nullptr,error);
        if (i % nclusterprocs) {
          iproc = me + (i % nclusterprocs);
          MPI_Send(&n,1,MPI_INT,iproc,0,world);
          MPI_Recv(&tmp,0,MPI_INT,iproc,0,world,MPI_STATUS_IGNORE);
          MPI_Rsend(buf,n,MPI_DOUBLE,iproc,0,world);
        }
      } else if (i % nclusterprocs == me - fileproc) {
        MPI_Recv(&n,1,MPI_INT,fileproc,0,world,MPI_STATUS_IGNORE);
        if (n > maxbuf) {
          maxbuf = n;
          memory->destroy(buf);
          memory->create(buf,maxbuf,"read_restart:buf");
        }
        MPI_Irecv(buf,n,MPI_DOUBLE,fileproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,fileproc,0,world);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
      }
      if (i % nclusterprocs == me - fileproc) {
        m = 0;
        while (m < n) m += avec->unpack_restart(&buf[m]);
      }
    }
    if (filereader && fp != nullptr) {
      fclose(fp);
      fp = nullptr;
    }
    MPI_Comm_free(&clustercomm);
  }
  delete[] file;
  memory->destroy(buf);
  if (multiproc || mpiioflag) {
    if (remapflag) {
      double **x = atom->x;
      imageint *image = atom->image;
      int nlocal = atom->nlocal;
      for (int i = 0; i < nlocal; i++)
        domain->remap(x[i],image[i]);
    }
    if (nextra)
      modify->add_fix(fmt::format("_read_restart all READ_RESTART {} {}",
                                  nextra,modify->nfix_restart_peratom));
    if (atom->map_style != Atom::MAP_NONE) {
      atom->map_init();
      atom->map_set();
    }
    if (domain->triclinic) domain->x2lamda(atom->nlocal);
    if (domain->triclinic) domain->lamda2x(atom->nlocal);
    if (nextra) {
      memory->destroy(atom->extra);
      memory->create(atom->extra,atom->nmax,nextra,"atom:extra");
      auto fix = dynamic_cast<FixReadRestart *>(modify->get_fix_by_id("_read_restart"));
      int *count = fix->count;
      double **extra = fix->extra;
      double **atom_extra = atom->extra;
      int nlocal = atom->nlocal;
      for (int i = 0; i < nlocal; i++)
        for (int j = 0; j < count[i]; j++)
          atom_extra[i][j] = extra[i][j];
      modify->delete_fix("_read_restart");
    }
  }
  bigint natoms;
  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (me == 0)
    utils::logmesg(lmp,"  {} atoms\n",natoms);
  if (natoms != atom->natoms)
    error->all(FLERR,"Did not assign all restart atoms correctly");
  atom->tag_check();
  if (atom->map_style != Atom::MAP_NONE) {
    atom->map_init();
    atom->map_set();
  }
  MPI_Barrier(world);
  if (comm->me == 0)
    utils::logmesg(lmp,"  read_restart CPU = {:.3f} seconds\n",platform::walltime()-time1);
  delete mpiio;
}
std::string ReadRestart::file_search(const std::string &inpfile)
{
  auto dirname = platform::path_dirname(inpfile);
  auto filename = platform::path_basename(inpfile);
  auto pattern = filename;
  auto loc = pattern.find('%');
  if (loc != std::string::npos) pattern.replace(loc,1,"base");
  bigint maxnum = -1;
  loc = pattern.find('*');
  if (loc != std::string::npos) {
    if (loc > 256)
      error->one(FLERR, "Filename part before '*' is too long to find restart with largest step");
    pattern.replace(loc,1,"\\d+");
    if (!platform::path_is_directory(dirname))
      error->one(FLERR,"Cannot open directory {} to search for restart file: {}",dirname);
    for (const auto &candidate : platform::list_directory(dirname)) {
      if (utils::strmatch(candidate,pattern)) {
        bigint num = ATOBIGINT(utils::strfind(candidate.substr(loc),"\\d+").c_str());
        if (num > maxnum) maxnum = num;
      }
    }
    if (maxnum < 0) error->one(FLERR,"Found no restart file matching pattern");
    filename.replace(filename.find('*'),1,std::to_string(maxnum));
  }
  return platform::path_join(dirname,filename);
}
void ReadRestart::header()
{
  int xperiodic(-1),yperiodic(-1),zperiodic(-1);
  int flag = read_int();
  while (flag >= 0) {
    if (flag == VERSION) {
      char *version = read_string();
      lmp->restart_ver = utils::date2num(version);
      if (me == 0)
        utils::logmesg(lmp,"  restart file = {}, LAMMPS = {}\n", version, lmp->version);
      delete[] version;
      if (revision > FORMAT_REVISION)
        error->all(FLERR,"Restart file format revision incompatible with current LAMMPS version");
      if ((me == 0) && (revision < FORMAT_REVISION))
        error->warning(FLERR,"Old restart file format revision. Switching to compatibility mode.");
    } else if (flag == SMALLINT) {
      int size = read_int();
      if (size != sizeof(smallint))
        error->all(FLERR,"Smallint setting in lmptype.h is not compatible");
    } else if (flag == IMAGEINT) {
      int size = read_int();
      if (size != sizeof(imageint))
        error->all(FLERR,"Imageint setting in lmptype.h is not compatible");
    } else if (flag == TAGINT) {
      int size = read_int();
      if (size != sizeof(tagint))
        error->all(FLERR,"Tagint setting in lmptype.h is not compatible");
    } else if (flag == BIGINT) {
      int size = read_int();
      if (size != sizeof(bigint))
        error->all(FLERR,"Bigint setting in lmptype.h is not compatible");
    } else if (flag == UNITS) {
      char *style = read_string();
      if (strcmp(style,update->unit_style) != 0) update->set_units(style);
      delete[] style;
    } else if (flag == NTIMESTEP) {
      update->ntimestep = read_bigint();
    } else if (flag == DIMENSION) {
      int dimension = read_int();
      domain->dimension = dimension;
      if (domain->dimension == 2 && domain->zperiodic == 0)
        error->all(FLERR, "Cannot run 2d simulation with non-periodic Z dimension");
    } else if (flag == NPROCS) {
      nprocs_file = read_int();
      if (nprocs_file != comm->nprocs && me == 0)
        error->warning(FLERR,"Restart file used different # of processors: {} vs. {}",
                       nprocs_file,comm->nprocs);
    } else if (flag == PROCGRID) {
      int procgrid[3];
      read_int();
      read_int_vec(3,procgrid);
      flag = 0;
      if (comm->user_procgrid[0] != 0 &&
          procgrid[0] != comm->user_procgrid[0]) flag = 1;
      if (comm->user_procgrid[1] != 0 &&
          procgrid[1] != comm->user_procgrid[1]) flag = 1;
      if (comm->user_procgrid[2] != 0 &&
          procgrid[2] != comm->user_procgrid[2]) flag = 1;
      if (flag && me == 0)
        error->warning(FLERR,"Restart file used different 3d processor grid");
    } else if (flag == NEWTON_PAIR) {
      int newton_pair_file = read_int();
      if (force->newton_pair != 1) {
        if (newton_pair_file != force->newton_pair && me == 0)
          error->warning(FLERR, "Restart file used different newton pair setting, "
                         "using input script value");
      }
    } else if (flag == XPERIODIC) {
      xperiodic = read_int();
    } else if (flag == YPERIODIC) {
      yperiodic = read_int();
    } else if (flag == ZPERIODIC) {
      zperiodic = read_int();
    } else if (flag == BOUNDARY) {
      int boundary[3][2];
      read_int();
      read_int_vec(6,&boundary[0][0]);
      if (domain->boundary[0][0] || domain->boundary[0][1] ||
          domain->boundary[1][0] || domain->boundary[1][1] ||
          domain->boundary[2][0] || domain->boundary[2][1]) {
        if (boundary[0][0] != domain->boundary[0][0] ||
            boundary[0][1] != domain->boundary[0][1] ||
            boundary[1][0] != domain->boundary[1][0] ||
            boundary[1][1] != domain->boundary[1][1] ||
            boundary[2][0] != domain->boundary[2][0] ||
            boundary[2][1] != domain->boundary[2][1]) {
          if (me == 0)
            error->warning(FLERR, "Restart file used different boundary settings, "
                           "using restart file values");
        }
      }
      domain->boundary[0][0] = boundary[0][0];
      domain->boundary[0][1] = boundary[0][1];
      domain->boundary[1][0] = boundary[1][0];
      domain->boundary[1][1] = boundary[1][1];
      domain->boundary[2][0] = boundary[2][0];
      domain->boundary[2][1] = boundary[2][1];
      if (xperiodic < 0 || yperiodic < 0 || zperiodic < 0)
        error->all(FLERR,"Illegal or unset periodicity in restart");
      domain->periodicity[0] = domain->xperiodic = xperiodic;
      domain->periodicity[1] = domain->yperiodic = yperiodic;
      domain->periodicity[2] = domain->zperiodic = zperiodic;
      domain->nonperiodic = 0;
      if (xperiodic == 0 || yperiodic == 0 || zperiodic == 0) {
        domain->nonperiodic = 1;
        if (boundary[0][0] >= 2 || boundary[0][1] >= 2 ||
            boundary[1][0] >= 2 || boundary[1][1] >= 2 ||
            boundary[2][0] >= 2 || boundary[2][1] >= 2)
          domain->nonperiodic = 2;
      }
    } else if (flag == BOUNDMIN) {
      double minbound[6];
      read_int();
      read_double_vec(6,minbound);
      domain->minxlo = minbound[0]; domain->minxhi = minbound[1];
      domain->minylo = minbound[2]; domain->minyhi = minbound[3];
      domain->minzlo = minbound[4]; domain->minzhi = minbound[5];
    } else if (flag == ATOM_STYLE) {
      char *style = read_string();
      int nargcopy = read_int();
      auto argcopy = new char*[nargcopy];
      for (int i = 0; i < nargcopy; i++)
        argcopy[i] = read_string();
      atom->create_avec(style,nargcopy,argcopy,1);
      if (comm->me ==0)
        utils::logmesg(lmp,"  restoring atom style {} from restart\n",atom->atom_style);
      for (int i = 0; i < nargcopy; i++) delete[] argcopy[i];
      delete[] argcopy;
      delete[] style;
    } else if (flag == NATOMS) {
      atom->natoms = read_bigint();
    } else if (flag == NTYPES) {
      atom->ntypes = read_int();
    } else if (flag == TRICLINIC) {
      domain->triclinic = read_int();
    } else if (flag == BOXLO) {
      read_int();
      read_double_vec(3,domain->boxlo);
    } else if (flag == BOXHI) {
      read_int();
      read_double_vec(3,domain->boxhi);
    } else if (flag == XY) {
      domain->xy = read_double();
    } else if (flag == XZ) {
      domain->xz = read_double();
    } else if (flag == YZ) {
      domain->yz = read_double();
    } else if (flag == TIMESTEP) {
      update->dt = read_double();
    } else if (flag == ATOM_ID) {
      atom->tag_enable = read_int();
    } else if (flag == ATOM_MAP_STYLE) {
      atom->map_style = read_int();
    } else if (flag == ATOM_MAP_USER) {
      atom->map_user = read_int();
    } else if (flag == ATOM_SORTFREQ) {
      atom->sortfreq = read_int();
    } else if (flag == ATOM_SORTBIN) {
      atom->userbinsize = read_double();
    } else if (flag == COMM_MODE) {
      comm->mode = read_int();
    } else if (flag == COMM_CUTOFF) {
      comm->cutghostuser = read_double();
    } else if (flag == COMM_VEL) {
      comm->ghost_velocity = read_int();
    } else if (flag == ATIMESTEP) {
      update->atimestep = read_bigint();
    } else if (flag == ATIME) {
      update->atime = read_double();
    } else error->all(FLERR,"Invalid flag in header section of restart file");
    flag = read_int();
  }
}
void ReadRestart::type_arrays()
{
  int flag = read_int();
  while (flag >= 0) {
    if (flag == MASS) {
      read_int();
      auto mass = new double[atom->ntypes+1];
      read_double_vec(atom->ntypes,&mass[1]);
      atom->set_mass(mass);
      delete[] mass;
    } else error->all(FLERR,
                      "Invalid flag in type arrays section of restart file");
    flag = read_int();
  }
}
void ReadRestart::file_layout()
{
  int flag = read_int();
  while (flag >= 0) {
    if (flag == MULTIPROC) {
      multiproc_file = read_int();
      if (multiproc == 0 && multiproc_file)
        error->all(FLERR,"Restart file is not a multi-proc file");
      if (multiproc && multiproc_file == 0)
        error->all(FLERR,"Restart file is a multi-proc file");
    } else if (flag == MPIIO) {
      int mpiioflag_file = read_int();
      if (mpiioflag == 0 && mpiioflag_file)
        error->all(FLERR,"Restart file is a MPI-IO file");
      if (mpiioflag && mpiioflag_file == 0)
        error->all(FLERR,"Restart file is not a MPI-IO file");
      if (mpiioflag) {
        bigint *nproc_chunk_offsets;
        memory->create(nproc_chunk_offsets,nprocs,
                       "write_restart:nproc_chunk_offsets");
        bigint *nproc_chunk_sizes;
        memory->create(nproc_chunk_sizes,nprocs,
                       "write_restart:nproc_chunk_sizes");
        if (me == 0) {
          int ndx;
          int *all_written_send_sizes;
          memory->create(all_written_send_sizes,nprocs_file,
                         "write_restart:all_written_send_sizes");
          int *nproc_chunk_number;
          memory->create(nproc_chunk_number,nprocs,
                         "write_restart:nproc_chunk_number");
          utils::sfread(FLERR,all_written_send_sizes,sizeof(int),nprocs_file,fp,nullptr,error);
          if ((nprocs != nprocs_file) && !(atom->nextra_store)) {
            atom->nlocal = 1;
            int perAtomSize = atom->avec->size_restart();
            atom->nlocal = 0;
            bigint total_size = 0;
            for (int i = 0; i < nprocs_file; ++i) {
              total_size += all_written_send_sizes[i];
            }
            bigint total_ct = total_size / perAtomSize;
            bigint base_ct = total_ct / nprocs;
            bigint leftover_ct = total_ct - (base_ct * nprocs);
            bigint current_ByteOffset = 0;
            base_ct += 1;
            bigint base_ByteOffset = base_ct * (perAtomSize * sizeof(double));
            for (ndx = 0; ndx < leftover_ct; ++ndx) {
              nproc_chunk_offsets[ndx] = current_ByteOffset;
              nproc_chunk_sizes[ndx] = base_ct * perAtomSize;
              current_ByteOffset += base_ByteOffset;
            }
            base_ct -= 1;
            base_ByteOffset -= (perAtomSize * sizeof(double));
            for (; ndx < nprocs; ++ndx) {
              nproc_chunk_offsets[ndx] = current_ByteOffset;
              nproc_chunk_sizes[ndx] = base_ct * perAtomSize;
              current_ByteOffset += base_ByteOffset;
            }
          } else {
            int init_chunk_number = nprocs_file/nprocs;
            int num_extra_chunks = nprocs_file - (nprocs*init_chunk_number);
            for (int i = 0; i < nprocs; i++) {
              if (i < num_extra_chunks)
                nproc_chunk_number[i] = init_chunk_number+1;
              else
                nproc_chunk_number[i] = init_chunk_number;
            }
            int all_written_send_sizes_index = 0;
            bigint current_offset = 0;
            for (int i=0;i<nprocs;i++) {
              nproc_chunk_offsets[i] = current_offset;
              nproc_chunk_sizes[i] = 0;
              for (int j=0;j<nproc_chunk_number[i];j++) {
                nproc_chunk_sizes[i] +=
                  all_written_send_sizes[all_written_send_sizes_index];
                current_offset +=
                  (all_written_send_sizes[all_written_send_sizes_index] *
                   sizeof(double));
                all_written_send_sizes_index++;
              }
            }
          }
          memory->destroy(all_written_send_sizes);
          memory->destroy(nproc_chunk_number);
        }
        MPI_Scatter(nproc_chunk_sizes, 1, MPI_LMP_BIGINT,
                    &assignedChunkSize , 1, MPI_LMP_BIGINT, 0,world);
        MPI_Scatter(nproc_chunk_offsets, 1, MPI_LMP_BIGINT,
                    &assignedChunkOffset , 1, MPI_LMP_BIGINT, 0,world);
        memory->destroy(nproc_chunk_sizes);
        memory->destroy(nproc_chunk_offsets);
      }
    }
    flag = read_int();
  }
  if (mpiioflag) {
    if (me == 0) headerOffset = platform::ftell(fp);
    MPI_Bcast(&headerOffset,1,MPI_LMP_BIGINT,0,world);
  }
}
void ReadRestart::magic_string()
{
  int n = strlen(MAGIC_STRING) + 1;
  auto str = new char[n];
  int count;
  if (me == 0) count = fread(str,sizeof(char),n,fp);
  MPI_Bcast(&count,1,MPI_INT,0,world);
  if (count < n)
    error->all(FLERR,"Invalid LAMMPS restart file");
  MPI_Bcast(str,n,MPI_CHAR,0,world);
  if (strcmp(str,MAGIC_STRING) != 0)
    error->all(FLERR,"Invalid LAMMPS restart file");
  delete[] str;
}
void ReadRestart::endian()
{
  int endian = read_int();
  if (endian == ENDIAN) return;
  if (endian == ENDIANSWAP)
    error->all(FLERR,"Restart file byte ordering is swapped");
  else error->all(FLERR,"Restart file byte ordering is not recognized");
}
void ReadRestart::format_revision()
{
  revision = read_int();
}
void ReadRestart::check_eof_magic()
{
  if (revision < 1) return;
  int n = strlen(MAGIC_STRING) + 1;
  auto str = new char[n];
  if (me == 0) {
    bigint curpos = platform::ftell(fp);
    platform::fseek(fp,platform::END_OF_FILE);
    bigint offset = platform::ftell(fp) - n;
    platform::fseek(fp,offset);
    utils::sfread(FLERR,str,sizeof(char),n,fp,nullptr,error);
    platform::fseek(fp,curpos);
  }
  MPI_Bcast(str,n,MPI_CHAR,0,world);
  if (strcmp(str,MAGIC_STRING) != 0)
    error->all(FLERR,"Incomplete or corrupted LAMMPS restart file");
  delete[] str;
}
int ReadRestart::read_int()
{
  int value;
  if ((me == 0) && (fread(&value,sizeof(int),1,fp) < 1))
    value = -1;
  MPI_Bcast(&value,1,MPI_INT,0,world);
  return value;
}
bigint ReadRestart::read_bigint()
{
  bigint value;
  if ((me == 0) && (fread(&value,sizeof(bigint),1,fp) < 1))
    value = -1;
  MPI_Bcast(&value,1,MPI_LMP_BIGINT,0,world);
  return value;
}
double ReadRestart::read_double()
{
  double value;
  if ((me == 0) && (fread(&value,sizeof(double),1,fp) < 1))
    value = 0.0;
  MPI_Bcast(&value,1,MPI_DOUBLE,0,world);
  return value;
}
char *ReadRestart::read_string()
{
  int n = read_int();
  if (n < 0) error->all(FLERR,"Illegal size string or corrupt restart");
  auto value = new char[n];
  if (me == 0) utils::sfread(FLERR,value,sizeof(char),n,fp,nullptr,error);
  MPI_Bcast(value,n,MPI_CHAR,0,world);
  return value;
}
void ReadRestart::read_int_vec(int n, int *vec)
{
  if (n < 0) error->all(FLERR,"Illegal size integer vector read requested");
  if (me == 0) utils::sfread(FLERR,vec,sizeof(int),n,fp,nullptr,error);
  MPI_Bcast(vec,n,MPI_INT,0,world);
}
void ReadRestart::read_double_vec(int n, double *vec)
{
  if (n < 0) error->all(FLERR,"Illegal size double vector read requested");
  if (me == 0) utils::sfread(FLERR,vec,sizeof(double),n,fp,nullptr,error);
  MPI_Bcast(vec,n,MPI_DOUBLE,0,world);
}
