#ifndef LMP_EXCEPTIONS_H
#define LMP_EXCEPTIONS_H 
#include <exception>
#include <mpi.h>
#include <string>
namespace LAMMPS_NS {
class LAMMPSException : public std::exception {
 public:
  std::string message;
  LAMMPSException(const std::string &msg) : message(msg) {}
  const char *what() const noexcept override { return message.c_str(); }
};
class LAMMPSAbortException : public LAMMPSException {
 public:
  MPI_Comm universe;
  LAMMPSAbortException(const std::string &msg, MPI_Comm _universe) :
      LAMMPSException(msg), universe(_universe)
  {
  }
};
enum ErrorType { ERROR_NONE = 0, ERROR_NORMAL = 1, ERROR_ABORT = 2 };
}
#endif
