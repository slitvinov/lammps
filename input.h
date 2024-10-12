#ifndef LMP_INPUT_H
#define LMP_INPUT_H 
#include "pointers.h"
#include <map>
namespace LAMMPS_NS {
class Command;
class Input : protected Pointers {
  friend class Info;
  friend class Error;
  friend class Deprecated;
  friend class SimpleCommandsTest_Echo_Test;
 public:
  int narg;
  char **arg;
  class Variable *variable;
  Input(class LAMMPS *, int, char **);
  ~Input() override;
  void file();
  void substitute(char *&, char *&, int &, int &, int);
  void write_echo(const std::string &);
 protected:
  char *command;
  int echo_screen;
  int echo_log;
 private:
  int me;
  int maxarg;
  char *line, *copy, *work;
  int maxline, maxcopy, maxwork;
  int nfile, maxfile;
  int label_active;
  char *labelstr;
  int jump_skip;
  bool utf8_warn;
  FILE **infiles;
 public:
  typedef Command *(*CommandCreator)(LAMMPS *);
  typedef std::map<std::string, CommandCreator> CommandCreatorMap;
  CommandCreatorMap *command_map;
 private:
  void parse();
  char *nextword(char *, char **);
  int numtriple(char *);
  void reallocate(char *&, int &, int);
  int execute_command();
  int meta(const std::string &);
  void log();
  void partition();
  void print();
  void quit();
  void atom_modify();
  void atom_style();
  void boundary();
  void comm_modify();
  void comm_style();
  void compute();
  void compute_modify();
  void dielectric();
  void dimension();
  void fix();
  void fix_modify();
  void group_command();
  void lattice();
  void mass();
  void neigh_modify();
  void neighbor_command();
  void newton();
  void pair_coeff();
  void pair_modify();
  void pair_style();
  void pair_write();
  void processors();
  void region();
  void run_style();
  void timestep();
  void uncompute();
  void unfix();
  void units();
};
}
#endif
