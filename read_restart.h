#ifdef COMMAND_CLASS
CommandStyle(read_restart,ReadRestart);
#else
#ifndef LMP_READ_RESTART_H
#define LMP_READ_RESTART_H 
#include "command.h"
namespace LAMMPS_NS {
class ReadRestart : public Command {
 public:
  ReadRestart(class LAMMPS *);
  void command(int, char **) override;
 private:
  int me, nprocs;
  FILE *fp;
  int multiproc;
  int multiproc_file;
  int nprocs_file;
  int revision;
  int mpiioflag;
  class RestartMPIIO *mpiio;
  bigint assignedChunkSize;
  MPI_Offset assignedChunkOffset, headerOffset;
  std::string file_search(const std::string &);
  void header();
  void type_arrays();
  void magic_string();
  void endian();
  void format_revision();
  void check_eof_magic();
  void file_layout();
  int read_int();
  bigint read_bigint();
  double read_double();
  char *read_string();
  void read_int_vec(int, int *);
  void read_double_vec(int, double *);
};
}
#endif
#endif
