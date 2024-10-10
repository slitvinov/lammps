/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

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
  int narg;                    // # of command args
  char **arg;                  // parsed args for command
  class Variable *variable;    // defined variables

  Input(class LAMMPS *, int, char **);
  ~Input() override;
  void file();                                             // process all input
  void file(const char *);                                 // process an input script
  char *one(const std::string &);                          // process a single command
  void substitute(char *&, char *&, int &, int &, int);    // substitute for variables in a string
  void write_echo(const std::string &);                    // send text to active echo file pointers

 protected:
  char *command;      // ptr to current command
  int echo_screen;    // 0 = no, 1 = yes
  int echo_log;       // 0 = no, 1 = yes

 private:
  int me;                           // proc ID
  int maxarg;                       // max # of args in arg
  char *line, *copy, *work;         // input line & copy and work string
  int maxline, maxcopy, maxwork;    // max lengths of char strings
  int nfile, maxfile;               // current # and max # of open input files
  int label_active;                 // 0 = no label, 1 = looking for label
  char *labelstr;                   // label string being looked for
  int jump_skip;                    // 1 if skipping next jump, 0 otherwise
  bool utf8_warn;                   // true if need to warn about UTF-8 chars

  FILE **infiles;    // list of open input files

 public:
  typedef Command *(*CommandCreator)(LAMMPS *);
  typedef std::map<std::string, CommandCreator> CommandCreatorMap;
  CommandCreatorMap *command_map;

 private:
  void parse();                            // parse an input text line
  char *nextword(char *, char **);         // find next word in string with quotes
  int numtriple(char *);                   // count number of triple quotes
  void reallocate(char *&, int &, int);    // reallocate a char string
  int execute_command();                   // execute a single command

  int meta(const std::string &);    // process meta-commands

  void clear();    // input script commands
  void echo();
  void ifthenelse();
  void include();
  void jump();
  void label();
  void log();
  void next_command();
  void partition();
  void plugin();
  void print();
  void quit();
  void shell();
  void variable_command();

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
  void package();
  void pair_coeff();
  void pair_modify();
  void pair_style();
  void pair_write();
  void processors();
  void region();
  void reset_timestep();
  void run_style();
  void suffix();
  void timestep();
  void uncompute();
  void unfix();
  void units();
};
}    // namespace LAMMPS_NS
#endif
