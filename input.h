#ifndef LMP_INPUT_H
#define LMP_INPUT_H
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
  Input(class LAMMPS *, int, char **);
  void file();
  void substitute(char *&, char *&, int &, int &, int);

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
};
} // namespace LAMMPS_NS
#endif
