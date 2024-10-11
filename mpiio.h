#ifndef LMP_MPIIO_H
#define LMP_MPIIO_H 
#ifdef LMP_MPIIO
#if defined(MPI_STUBS)
#error "The MPIIO package cannot be compiled in serial with MPI STUBS"
#endif
#include "restart_mpiio.h"
#else
namespace LAMMPS_NS {
class RestartMPIIO {
 public:
  int mpiio_exists;
  RestartMPIIO(class LAMMPS *) { mpiio_exists = 0; }
  ~RestartMPIIO() {}
  void openForRead(const char *) {}
  void openForWrite(const char *) {}
  void write(MPI_Offset, int, double *) {}
  void read(MPI_Offset, long, double *) {}
  void close() {}
};
}
#endif
#endif
