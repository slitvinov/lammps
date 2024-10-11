#ifdef COMMAND_CLASS
CommandStyle(balance,Balance);
#else
#ifndef LMP_BALANCE_H
#define LMP_BALANCE_H 
#include "command.h"
namespace LAMMPS_NS {
class Balance : public Command {
 public:
  class RCB *rcb;
  int sortflag;
  int outflag;
  Balance(class LAMMPS *);
  ~Balance() override;
  void command(int, char **) override;
  void options(int, int, char **, int);
  double imbalance_factor(double &);
  void shift_setup(char *, int, double);
  int shift();
  int *bisection();
  void dumpout(bigint);
  static constexpr int BSTR_SIZE = 3;
 private:
  int me, nprocs;
  double thresh;
  int style;
  int xflag, yflag, zflag;
  double *user_xsplit, *user_ysplit, *user_zsplit;
  int oldrcb;
  int nitermax;
  double stopthresh;
  char bstr[BSTR_SIZE + 1];
  int shift_allocate;
  int ndim;
  int *bdim;
  double *onecost;
  double *allcost;
  double *sum;
  double *target;
  double *lo, *hi;
  double *losum, *hisum;
  int rho;
  double *proccost;
  double *allproccost;
  int nimbalance;
  class Imbalance **imbalances;
  FILE *fp;
  int firststep;
  double imbalance_splits();
  void shift_setup_static(char *);
  void tally(int, int, double *);
  int adjust(int, double *);
#ifdef BALANCE_DEBUG
  void debug_shift_output(int, int, int, double *);
#endif
};
}
#endif
#endif
