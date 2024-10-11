#include "variable.h"
#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "group.h"
#include "input.h"
#include "label_map.h"
#include "library.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "random_mars.h"
#include "region.h"
#include "tokenizer.h"
#include "universe.h"
#include "update.h"
#include <cctype>
#include <cmath>
#include <cstring>
#include <unordered_map>
using namespace LAMMPS_NS;
using namespace MathConst;
#define VARDELTA 4
#define MAXLEVEL 4
#define MAXLINE 256
#define CHUNK 1024
#define MAXFUNCARG 6
#define MYROUND(a) (( (a)-floor(a) ) >= .5) ? ceil(a) : floor(a)
enum{ARG,OP};
enum{DONE,ADD,SUBTRACT,MULTIPLY,DIVIDE,CARAT,MODULO,UNARY,
     NOT,EQ,NE,LT,LE,GT,GE,AND,OR,XOR,
     SQRT,EXP,LN,LOG,ABS,SIN,COS,TAN,ASIN,ACOS,ATAN,ATAN2,
     RANDOM,NORMAL,CEIL,FLOOR,ROUND,RAMP,STAGGER,LOGFREQ,LOGFREQ2,
     LOGFREQ3,STRIDE,STRIDE2,VDISPLACE,SWIGGLE,CWIGGLE,GMASK,RMASK,
     GRMASK,IS_ACTIVE,IS_DEFINED,IS_AVAILABLE,IS_FILE,EXTRACT_SETTING,
     VALUE,ATOMARRAY,TYPEARRAY,INTARRAY,BIGINTARRAY,VECTORARRAY};
enum{SUM,XMIN,XMAX,AVE,TRAP,SLOPE};
static constexpr double BIG = 1.0e20;
#if defined(LAMMPS_SMALLBIG) || defined(LAMMPS_BIGBIG)
static constexpr double MAXBIGINT_DOUBLE = (double) (MAXBIGINT-512);
#else
static constexpr double MAXBIGINT_DOUBLE = (double) MAXBIGINT;
#endif
static std::unordered_map<std::string, double> constants = {
  {"PI", MY_PI },
  {"version", -1 },
  {"yes", 1 },
  {"no", 0 },
  {"on", 1 },
  {"off", 0 },
  {"true", 1 },
  {"false", 0 }
};
Variable::Variable(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  nvar = maxvar = 0;
  names = nullptr;
  style = nullptr;
  num = nullptr;
  which = nullptr;
  pad = nullptr;
  data = nullptr;
  dvalue = nullptr;
  vecs = nullptr;
  eval_in_progress = nullptr;
  randomequal = nullptr;
  randomatom = nullptr;
  constants["version"] = lmp->num_ver;
  precedence[DONE] = 0;
  precedence[OR] = precedence[XOR] = 1;
  precedence[AND] = 2;
  precedence[EQ] = precedence[NE] = 3;
  precedence[LT] = precedence[LE] = precedence[GT] = precedence[GE] = 4;
  precedence[ADD] = precedence[SUBTRACT] = 5;
  precedence[MULTIPLY] = precedence[DIVIDE] = precedence[MODULO] = 6;
  precedence[CARAT] = 7;
  precedence[UNARY] = precedence[NOT] = 8;
}
Variable::~Variable()
{
  for (int i = 0; i < nvar; i++) {
    delete[] names[i];
    if (style[i] == LOOP || style[i] == ULOOP) delete[] data[i][0];
    else for (int j = 0; j < num[i]; j++) delete[] data[i][j];
    delete[] data[i];
    if (style[i] == VECTOR) memory->destroy(vecs[i].values);
  }
  memory->sfree(names);
  memory->destroy(style);
  memory->destroy(num);
  memory->destroy(which);
  memory->destroy(pad);
  memory->sfree(data);
  memory->sfree(dvalue);
  memory->sfree(vecs);
  memory->destroy(eval_in_progress);
  delete randomequal;
  delete randomatom;
}
void Variable::set(int narg, char **arg)
{
  if (narg < 2) utils::missing_cmd_args(FLERR, "variable", error);
  int replaceflag = 0;
  if (strcmp(arg[1],"equal") == 0) {
    if (narg != 3) error->all(FLERR,"Illegal variable command: expected 3 arguments but found {}", narg);
    int ivar = find(arg[0]);
    if (ivar >= 0) {
      if (style[ivar] != EQUAL)
        error->all(FLERR,"Cannot redefine variable as a different style");
      delete[] data[ivar][0];
      data[ivar][0] = utils::strdup(arg[2]);
      replaceflag = 1;
    } else {
      if (nvar == maxvar) grow();
      style[nvar] = EQUAL;
      num[nvar] = 2;
      which[nvar] = 0;
      pad[nvar] = 0;
      data[nvar] = new char*[num[nvar]];
      data[nvar][0] = utils::strdup(arg[2]);
      data[nvar][1] = new char[VALUELENGTH];
      strcpy(data[nvar][1],"(undefined)");
    }
  } else error->all(FLERR,"Unknown variable keyword: {}", arg[1]);
  if (replaceflag) return;
  if (!utils::is_id(arg[0]))
    error->all(FLERR,"Variable name '{}' must have only letters, numbers, or underscores",arg[0]);
  names[nvar] = utils::strdup(arg[0]);
  nvar++;
}
int Variable::find(const char *name)
{
  if (name == nullptr) return -1;
  for (int i = 0; i < nvar; i++)
    if (strcmp(name,names[i]) == 0) return i;
  return -1;
}
int Variable::equalstyle(int ivar)
{
  if (style[ivar] == EQUAL || style[ivar] == TIMER ||
      style[ivar] == INTERNAL) return 1;
  return 0;
}
int Variable::atomstyle(int ivar)
{
  if (style[ivar] == ATOM ) return 1;
  return 0;
}
int Variable::vectorstyle(int ivar)
{
  if (style[ivar] == VECTOR) return 1;
  return 0;
}
int Variable::internalstyle(int ivar)
{
  if (style[ivar] == INTERNAL) return 1;
  return 0;
}
char *Variable::retrieve(const char *name)
{
  int ivar = find(name);
  if (ivar < 0) return nullptr;
  if (which[ivar] >= num[ivar]) return nullptr;
  if (eval_in_progress[ivar])
    print_var_error(FLERR,"has a circular dependency",ivar);
  eval_in_progress[ivar] = 1;
  char *str = nullptr;
  if (style[ivar] == INDEX || style[ivar] == WORLD ||
      style[ivar] == UNIVERSE || style[ivar] == STRING ||
      style[ivar] == SCALARFILE) {
    str = data[ivar][which[ivar]];
  } else if (style[ivar] == LOOP || style[ivar] == ULOOP) {
    std::string result;
    if (pad[ivar] == 0) result = std::to_string(which[ivar]+1);
    else result = fmt::format("{:0>{}d}",which[ivar]+1, pad[ivar]);
    delete[] data[ivar][0];
    str = data[ivar][0] = utils::strdup(result);
  } else if (style[ivar] == EQUAL) {
    double answer = evaluate(data[ivar][0],nullptr,ivar);
    delete[] data[ivar][1];
    data[ivar][1] = utils::strdup(fmt::format("{:.15g}",answer));
    str = data[ivar][1];
  } else if (style[ivar] == FORMAT) {
    int jvar = find(data[ivar][0]);
    if (jvar < 0)
      error->all(FLERR, "Variable {}: format variable {} does not exist", names[ivar],data[ivar][0]);
    if (!equalstyle(jvar))
      error->all(FLERR, "Variable {}: format variable {} has incompatible style",
                 names[ivar],data[ivar][0]);
    double answer = compute_equal(jvar);
    sprintf(data[ivar][2],data[ivar][1],answer);
    str = data[ivar][2];
  } else if (style[ivar] == GETENV) {
    const char *result = getenv(data[ivar][0]);
    if (result == nullptr) result = (const char *) "";
    delete[] data[ivar][1];
    str = data[ivar][1] = utils::strdup(result);
  } else if (style[ivar] == TIMER || style[ivar] == INTERNAL) {
    delete[] data[ivar][0];
    data[ivar][0] = utils::strdup(fmt::format("{:.15g}",dvalue[ivar]));
    str = data[ivar][0];
  } else if (style[ivar] == ATOM ||
             style[ivar] == VECTOR) return nullptr;
  eval_in_progress[ivar] = 0;
  return str;
}
double Variable::compute_equal(int ivar)
{
  if (eval_in_progress[ivar])
    print_var_error(FLERR,"has a circular dependency",ivar);
  eval_in_progress[ivar] = 1;
  double value = 0.0;
  if (style[ivar] == EQUAL) value = evaluate(data[ivar][0],nullptr,ivar);
  else if (style[ivar] == TIMER) value = dvalue[ivar];
  else if (style[ivar] == INTERNAL) value = dvalue[ivar];
  eval_in_progress[ivar] = 0;
  return value;
}
double Variable::compute_equal(const std::string &str)
{
  char *ptr = utils::strdup(str);
  double val = evaluate(ptr,nullptr,-1);
  delete[] ptr;
  return val;
}
void Variable::compute_atom(int ivar, int igroup, double *result, int stride, int sumflag)
{
  Tree *tree = nullptr;
  double *vstore;
  if (eval_in_progress[ivar])
    print_var_error(FLERR,"has a circular dependency",ivar);
  eval_in_progress[ivar] = 1;
    treetype = ATOM;
    evaluate(data[ivar][0],&tree,ivar);
    collapse_tree(tree);
  if (result == nullptr) {
    if (style[ivar] == ATOM) free_tree(tree);
    eval_in_progress[ivar] = 0;
    return;
  }
  int groupbit = group->bitmask[igroup];
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (style[ivar] == ATOM) {
    if (sumflag == 0) {
      int m = 0;
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) result[m] = eval_tree(tree,i);
        else result[m] = 0.0;
        m += stride;
      }
    } else {
      int m = 0;
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) result[m] += eval_tree(tree,i);
        m += stride;
      }
    }
  } else {
    if (sumflag == 0) {
      int m = 0;
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) result[m] = vstore[i];
        else result[m] = 0.0;
        m += stride;
      }
    } else {
      int m = 0;
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) result[m] += vstore[i];
        m += stride;
      }
    }
  }
  if (style[ivar] == ATOM) free_tree(tree);
  eval_in_progress[ivar] = 0;
}
int Variable::compute_vector(int ivar, double **result)
{
  Tree *tree = nullptr;
  if (vecs[ivar].currentstep == update->ntimestep) {
    *result = vecs[ivar].values;
    return vecs[ivar].n;
  }
  if (eval_in_progress[ivar])
    print_var_error(FLERR,"has a circular dependency",ivar);
  eval_in_progress[ivar] = 1;
  treetype = VECTOR;
  evaluate(data[ivar][0],&tree,ivar);
  collapse_tree(tree);
  int nlen = size_tree_vector(tree);
  if (nlen == 0)
    print_var_error(FLERR,"Vector-style variable has zero length",ivar);
  if (nlen < 0)
    print_var_error(FLERR,"Inconsistent lengths in vector-style variable",ivar);
  if (nlen > vecs[ivar].nmax) {
    memory->destroy(vecs[ivar].values);
    vecs[ivar].nmax = nlen;
    memory->create(vecs[ivar].values,vecs[ivar].nmax,"variable:values");
  }
  vecs[ivar].n = nlen;
  vecs[ivar].currentstep = update->ntimestep;
  double *vec = vecs[ivar].values;
  for (int i = 0; i < nlen; i++)
    vec[i] = eval_tree(tree,i);
  free_tree(tree);
  eval_in_progress[ivar] = 0;
  *result = vec;
  return nlen;
}
void Variable::internal_set(int ivar, double value)
{
  dvalue[ivar] = value;
}
void Variable::remove(int n)
{
  delete[] names[n];
  if (style[n] == LOOP || style[n] == ULOOP) delete[] data[n][0];
  else for (int i = 0; i < num[n]; i++) delete[] data[n][i];
  delete[] data[n];
  for (int i = n+1; i < nvar; i++) {
    names[i-1] = names[i];
    style[i-1] = style[i];
    num[i-1] = num[i];
    which[i-1] = which[i];
    pad[i-1] = pad[i];
    data[i-1] = data[i];
    dvalue[i-1] = dvalue[i];
  }
  nvar--;
  data[nvar] = nullptr;
  names[nvar] = nullptr;
}
void Variable::grow()
{
  int old = maxvar;
  maxvar += VARDELTA;
  names = (char **) memory->srealloc(names,maxvar*sizeof(char *),"var:names");
  memory->grow(style,maxvar,"var:style");
  memory->grow(num,maxvar,"var:num");
  memory->grow(which,maxvar,"var:which");
  memory->grow(pad,maxvar,"var:pad");
  data = (char ***) memory->srealloc(data,maxvar*sizeof(char **),"var:data");
  memory->grow(dvalue,maxvar,"var:dvalue");
  vecs = (VecVar *) memory->srealloc(vecs,maxvar*sizeof(VecVar),"var:vecvar");
  for (int i = old; i < maxvar; i++) {
    vecs[i].nmax = 0;
    vecs[i].currentstep = -1;
    vecs[i].values = nullptr;
  }
  memory->grow(eval_in_progress,maxvar,"var:eval_in_progress");
  for (int i = 0; i < maxvar; i++) eval_in_progress[i] = 0;
}
void Variable::copy(int narg, char **from, char **to)
{
  for (int i = 0; i < narg; i++)
    to[i] = utils::strdup(from[i]);
}
double Variable::evaluate(char *str, Tree **tree, int ivar)
{
  int op,opprevious;
  double value1,value2;
  char onechar;
  char *ptr;
  double argstack[MAXLEVEL];
  Tree *treestack[MAXLEVEL];
  int opstack[MAXLEVEL];
  int nargstack = 0;
  int ntreestack = 0;
  int nopstack = 0;
  int i = 0;
  int expect = ARG;
  if (str == nullptr)
    print_var_error(FLERR,"Invalid syntax in variable formula",ivar);
  while (true) {
    onechar = str[i];
    if (isspace(onechar)) i++;
    else if (onechar == '(') {
      if (expect == OP)
        print_var_error(FLERR,"Invalid syntax in variable formula",ivar);
      expect = OP;
      char *contents = nullptr;
      i = find_matching_paren(str,i,contents,ivar);
      i++;
      if (tree) {
        Tree *newtree = nullptr;
        evaluate(contents,&newtree,ivar);
        treestack[ntreestack++] = newtree;
      } else argstack[nargstack++] = evaluate(contents,nullptr,ivar);
      delete[] contents;
    } else if (isdigit(onechar) || onechar == '.') {
      if (expect == OP)
        print_var_error(FLERR,"Invalid syntax in variable formula",ivar);
      expect = OP;
      int istart = i;
      while (isdigit(str[i]) || str[i] == '.') i++;
      if (str[i] == 'e' || str[i] == 'E') {
        i++;
        if (str[i] == '+' || str[i] == '-') i++;
        while (isdigit(str[i])) i++;
      }
      int istop = i - 1;
      int n = istop - istart + 1;
      auto number = new char[n+1];
      strncpy(number,&str[istart],n);
      number[n] = '\0';
      if (tree) {
        auto newtree = new Tree();
        newtree->type = VALUE;
        newtree->value = atof(number);
        treestack[ntreestack++] = newtree;
      } else argstack[nargstack++] = atof(number);
      delete[] number;
    } else if (isalpha(onechar)) {
      if (expect == OP)
        print_var_error(FLERR,"Invalid syntax in variable formula",ivar);
      expect = OP;
      int istart = i;
      while (isalnum(str[i]) || str[i] == '_') i++;
      int istop = i-1;
      int n = istop - istart + 1;
      auto word = new char[n+1];
      strncpy(word,&str[istart],n);
      word[n] = '\0';
      if (utils::strmatch(word,"^[Cc]_")) {
        if (domain->box_exist == 0)
          print_var_error(FLERR,"Variable evaluation before simulation box is defined",ivar);
        int lowercase = 1;
        if (word[0] == 'C') lowercase = 0;
        Compute *compute = modify->get_compute_by_id(word+2);
        if (!compute)
          print_var_error(FLERR,fmt::format("Invalid compute ID '{}' in variable formula", word+2),ivar);
        int nbracket;
        tagint index1,index2;
        if (str[i] != '[') nbracket = 0;
        else {
          nbracket = 1;
          ptr = &str[i];
          index1 = int_between_brackets(ptr,1);
          i = ptr-str+1;
          if (str[i] == '[') {
            nbracket = 2;
            ptr = &str[i];
            index2 = int_between_brackets(ptr,1);
            i = ptr-str+1;
          }
        }
        if (nbracket == 0 && compute->scalar_flag && lowercase) {
          if (update->whichflag == 0) {
            if (compute->invoked_scalar != update->ntimestep)
              print_var_error(FLERR,"Compute used in variable between runs is not current",ivar);
          } else if (!(compute->invoked_flag & Compute::INVOKED_SCALAR)) {
            compute->compute_scalar();
            compute->invoked_flag |= Compute::INVOKED_SCALAR;
          }
          value1 = compute->scalar;
          if (tree) {
            auto newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;
        } else if (nbracket == 1 && compute->vector_flag && lowercase) {
          if (index1 > compute->size_vector &&
              compute->size_vector_variable == 0)
            print_var_error(FLERR,"Variable formula compute vector is accessed out-of-range",ivar,0);
          if (update->whichflag == 0) {
            if (compute->invoked_vector != update->ntimestep)
              print_var_error(FLERR,"Compute used in variable between runs is not current",ivar);
          } else if (!(compute->invoked_flag & Compute::INVOKED_VECTOR)) {
            compute->compute_vector();
            compute->invoked_flag |= Compute::INVOKED_VECTOR;
          }
          if (compute->size_vector_variable &&
              index1 > compute->size_vector) value1 = 0.0;
          else value1 = compute->vector[index1-1];
          if (tree) {
            auto newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;
        } else if (nbracket == 2 && compute->array_flag && lowercase) {
          if (index1 > compute->size_array_rows &&
              compute->size_array_rows_variable == 0)
            print_var_error(FLERR,"Variable formula compute array is accessed out-of-range",ivar,0);
          if (index2 > compute->size_array_cols)
            print_var_error(FLERR,"Variable formula compute array is accessed out-of-range",ivar,0);
          if (update->whichflag == 0) {
            if (compute->invoked_array != update->ntimestep)
              print_var_error(FLERR,"Compute used in variable between runs is not current",ivar);
          } else if (!(compute->invoked_flag & Compute::INVOKED_ARRAY)) {
            compute->compute_array();
            compute->invoked_flag |= Compute::INVOKED_ARRAY;
          }
          if (compute->size_array_rows_variable &&
              index1 > compute->size_array_rows) value1 = 0.0;
          else value1 = compute->array[index1-1][index2-1];
          if (tree) {
            auto newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;
        } else if (nbracket == 0 && compute->vector_flag) {
          if (tree == nullptr)
            print_var_error(FLERR,"Compute global vector in equal-style variable formula",ivar);
          if (treetype == ATOM)
            print_var_error(FLERR,"Compute global vector in atom-style variable formula",ivar);
          if (compute->size_vector == 0)
            print_var_error(FLERR,"Variable formula compute vector is zero length",ivar);
          if (update->whichflag == 0) {
            if (compute->invoked_vector != update->ntimestep)
              print_var_error(FLERR,"Compute used in variable between runs is not current",ivar);
          } else if (!(compute->invoked_flag & Compute::INVOKED_VECTOR)) {
            compute->compute_vector();
            compute->invoked_flag |= Compute::INVOKED_VECTOR;
          }
          auto newtree = new Tree();
          newtree->type = VECTORARRAY;
          newtree->array = compute->vector;
          newtree->nvector = compute->size_vector;
          newtree->nstride = 1;
          treestack[ntreestack++] = newtree;
        } else if (nbracket == 1 && compute->array_flag) {
          if (tree == nullptr)
            print_var_error(FLERR,"Compute global vector in equal-style variable formula",ivar);
          if (treetype == ATOM)
            print_var_error(FLERR,"Compute global vector in atom-style variable formula",ivar);
          if (compute->size_array_rows == 0)
            print_var_error(FLERR,"Variable formula compute array is zero length",ivar);
          if (update->whichflag == 0) {
            if (compute->invoked_array != update->ntimestep)
              print_var_error(FLERR,"Compute used in variable between runs is not current",ivar);
          } else if (!(compute->invoked_flag & Compute::INVOKED_ARRAY)) {
            compute->compute_array();
            compute->invoked_flag |= Compute::INVOKED_ARRAY;
          }
          auto newtree = new Tree();
          newtree->type = VECTORARRAY;
          newtree->array = &compute->array[0][index1-1];
          newtree->nvector = compute->size_array_rows;
          newtree->nstride = compute->size_array_cols;
          treestack[ntreestack++] = newtree;
        } else if (nbracket == 1 && compute->peratom_flag &&
                   compute->size_peratom_cols == 0) {
          if (update->whichflag == 0) {
            if (compute->invoked_peratom != update->ntimestep)
              print_var_error(FLERR,"Compute used in variable between runs is not current",ivar);
          } else if (!(compute->invoked_flag & Compute::INVOKED_PERATOM)) {
            compute->compute_peratom();
            compute->invoked_flag |= Compute::INVOKED_PERATOM;
          }
          peratom2global(1,nullptr,compute->vector_atom,1,index1,tree,
                         treestack,ntreestack,argstack,nargstack);
        } else if (nbracket == 2 && compute->peratom_flag &&
                   compute->size_peratom_cols > 0) {
          if (index2 > compute->size_peratom_cols)
            print_var_error(FLERR,"Variable formula compute array is accessed out-of-range",ivar,0);
          if (update->whichflag == 0) {
            if (compute->invoked_peratom != update->ntimestep)
              print_var_error(FLERR,"Compute used in variable between runs is not current",ivar);
          } else if (!(compute->invoked_flag & Compute::INVOKED_PERATOM)) {
            compute->compute_peratom();
            compute->invoked_flag |= Compute::INVOKED_PERATOM;
          }
          if (compute->array_atom)
            peratom2global(1,nullptr,&compute->array_atom[0][index2-1],compute->size_peratom_cols,index1,
                           tree,treestack,ntreestack,argstack,nargstack);
          else
            peratom2global(1,nullptr,nullptr,compute->size_peratom_cols,index1,
                           tree,treestack,ntreestack,argstack,nargstack);
        } else if (nbracket == 0 && compute->peratom_flag &&
                   compute->size_peratom_cols == 0) {
          if (tree == nullptr)
            print_var_error(FLERR,"Per-atom compute in equal-style variable formula",ivar);
          if (treetype == VECTOR)
            print_var_error(FLERR,"Per-atom compute in vector-style variable formula",ivar);
          if (update->whichflag == 0) {
            if (compute->invoked_peratom != update->ntimestep)
              print_var_error(FLERR,"Compute used in variable between runs is not current",ivar);
          } else if (!(compute->invoked_flag & Compute::INVOKED_PERATOM)) {
            compute->compute_peratom();
            compute->invoked_flag |= Compute::INVOKED_PERATOM;
          }
          auto newtree = new Tree();
          newtree->type = ATOMARRAY;
          newtree->array = compute->vector_atom;
          newtree->nstride = 1;
          treestack[ntreestack++] = newtree;
        } else if (nbracket == 1 && compute->peratom_flag &&
                   compute->size_peratom_cols > 0) {
          if (tree == nullptr)
            print_var_error(FLERR,"Per-atom compute in equal-style variable formula",ivar);
          if (treetype == VECTOR)
            print_var_error(FLERR,"Per-atom compute in vector-style variable formula",ivar);
          if (index1 > compute->size_peratom_cols)
            print_var_error(FLERR,"Variable formula compute array is accessed out-of-range",ivar,0);
          if (update->whichflag == 0) {
            if (compute->invoked_peratom != update->ntimestep)
              print_var_error(FLERR,"Compute used in variable between runs is not current",ivar);
          } else if (!(compute->invoked_flag & Compute::INVOKED_PERATOM)) {
            compute->compute_peratom();
            compute->invoked_flag |= Compute::INVOKED_PERATOM;
          }
          auto newtree = new Tree();
          newtree->type = ATOMARRAY;
          if (compute->array_atom)
            newtree->array = &compute->array_atom[0][index1-1];
          newtree->nstride = compute->size_peratom_cols;
          treestack[ntreestack++] = newtree;
        } else if (nbracket == 1 && compute->local_flag) {
          print_var_error(FLERR,"Cannot access local data via indexing",ivar);
        } else print_var_error(FLERR,"Mismatched compute in variable formula",ivar);
      } else if (utils::strmatch(word,"^[fF]_")) {
        if (domain->box_exist == 0)
          print_var_error(FLERR,"Variable evaluation before simulation box is defined",ivar);
        int lowercase = 1;
        if (word[0] == 'F') lowercase = 0;
        Fix *fix = modify->get_fix_by_id(word+2);
        if (!fix)
          print_var_error(FLERR,fmt::format("Invalid fix ID '{}' in variable formula",word+2),ivar);
        int nbracket;
        tagint index1,index2;
        if (str[i] != '[') nbracket = 0;
        else {
          nbracket = 1;
          ptr = &str[i];
          index1 = int_between_brackets(ptr,1);
          i = ptr-str+1;
          if (str[i] == '[') {
            nbracket = 2;
            ptr = &str[i];
            index2 = int_between_brackets(ptr,1);
            i = ptr-str+1;
          }
        }
        if (nbracket == 0 && fix->scalar_flag && lowercase) {
          if (update->whichflag > 0 && update->ntimestep % fix->global_freq)
            print_var_error(FLERR,"Fix in variable not computed at a compatible time",ivar);
          value1 = fix->compute_scalar();
          if (tree) {
            auto newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;
        } else if (nbracket == 1 && fix->vector_flag && lowercase) {
          if (index1 > fix->size_vector &&
              fix->size_vector_variable == 0)
            print_var_error(FLERR,"Variable formula fix vector is accessed out-of-range",ivar,0);
          if (update->whichflag > 0 && update->ntimestep % fix->global_freq)
            print_var_error(FLERR,"Fix in variable not computed at a compatible time",ivar);
          value1 = fix->compute_vector(index1-1);
          if (tree) {
            auto newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;
        } else if (nbracket == 2 && fix->array_flag && lowercase) {
          if (index1 > fix->size_array_rows &&
              fix->size_array_rows_variable == 0)
            print_var_error(FLERR,"Variable formula fix array is accessed out-of-range",ivar,0);
          if (index2 > fix->size_array_cols)
            print_var_error(FLERR,"Variable formula fix array is accessed out-of-range",ivar,0);
          if (update->whichflag > 0 && update->ntimestep % fix->global_freq)
            print_var_error(FLERR,"Fix in variable not computed at a compatible time",ivar);
          value1 = fix->compute_array(index1-1,index2-1);
          if (tree) {
            auto newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;
        } else if (nbracket == 0 && fix->vector_flag) {
          if (update->whichflag > 0 && update->ntimestep % fix->global_freq)
            print_var_error(FLERR,"Fix in variable not computed at compatible time",ivar);
          if (tree == nullptr)
            print_var_error(FLERR,"Fix global vector in equal-style variable formula",ivar);
          if (treetype == ATOM)
            print_var_error(FLERR,"Fix global vector in atom-style variable formula",ivar);
          if (fix->size_vector == 0)
            print_var_error(FLERR,"Variable formula fix vector is zero length",ivar);
          int nvec = fix->size_vector;
          double *vec;
          memory->create(vec,nvec,"variable:values");
          for (int m = 0; m < nvec; m++)
            vec[m] = fix->compute_vector(m);
          auto newtree = new Tree();
          newtree->type = VECTORARRAY;
          newtree->array = vec;
          newtree->nvector = nvec;
          newtree->nstride = 1;
          newtree->selfalloc = 1;
          treestack[ntreestack++] = newtree;
        } else if (nbracket == 1 && fix->array_flag) {
          if (update->whichflag > 0 && update->ntimestep % fix->global_freq)
            print_var_error(FLERR,"Fix in variable not computed at a compatible time",ivar);
          if (tree == nullptr)
            print_var_error(FLERR,"Fix global vector in equal-style variable formula",ivar);
          if (treetype == ATOM)
            print_var_error(FLERR,"Fix global vector in atom-style variable formula",ivar);
          if (fix->size_array_rows == 0)
            print_var_error(FLERR,"Variable formula fix array is zero length",ivar);
          int nvec = fix->size_array_rows;
          double *vec;
          memory->create(vec,nvec,"variable:values");
          for (int m = 0; m < nvec; m++)
            vec[m] = fix->compute_array(m,index1-1);
          auto newtree = new Tree();
          newtree->type = VECTORARRAY;
          newtree->array = vec;
          newtree->nvector = nvec;
          newtree->nstride = 1;
          newtree->selfalloc = 1;
          treestack[ntreestack++] = newtree;
        } else if (nbracket == 1 && fix->peratom_flag &&
                   fix->size_peratom_cols == 0) {
          if (update->whichflag > 0 &&
              update->ntimestep % fix->peratom_freq)
            print_var_error(FLERR,"Fix in variable not computed at a compatible time",ivar);
          peratom2global(1,nullptr,fix->vector_atom,1,index1,
                         tree,treestack,ntreestack,argstack,nargstack);
        } else if (nbracket == 2 && fix->peratom_flag &&
                   fix->size_peratom_cols > 0) {
          if (index2 > fix->size_peratom_cols)
            print_var_error(FLERR,"Variable formula fix array is accessed out-of-range",ivar,0);
          if (update->whichflag > 0 &&
              update->ntimestep % fix->peratom_freq)
            print_var_error(FLERR,"Fix in variable not computed at a compatible time",ivar);
          if (fix->array_atom)
            peratom2global(1,nullptr,&fix->array_atom[0][index2-1],fix->size_peratom_cols,index1,
                           tree,treestack,ntreestack,argstack,nargstack);
          else
            peratom2global(1,nullptr,nullptr,fix->size_peratom_cols,index1,
                           tree,treestack,ntreestack,argstack,nargstack);
        } else if (nbracket == 0 && fix->peratom_flag &&
                   fix->size_peratom_cols == 0) {
          if (tree == nullptr)
            print_var_error(FLERR,"Per-atom fix in equal-style variable formula",ivar);
          if (update->whichflag > 0 &&
              update->ntimestep % fix->peratom_freq)
            print_var_error(FLERR,"Fix in variable not computed at compatible time",ivar);
          auto newtree = new Tree();
          newtree->type = ATOMARRAY;
          newtree->array = fix->vector_atom;
          newtree->nstride = 1;
          treestack[ntreestack++] = newtree;
        } else if (nbracket == 1 && fix->peratom_flag &&
                   fix->size_peratom_cols > 0) {
          if (tree == nullptr)
            print_var_error(FLERR,"Per-atom fix in equal-style variable formula",ivar);
          if (index1 > fix->size_peratom_cols)
            print_var_error(FLERR,"Variable formula fix array is accessed out-of-range",ivar,0);
          if (update->whichflag > 0 &&
              update->ntimestep % fix->peratom_freq)
            print_var_error(FLERR,"Fix in variable not computed at compatible time",ivar);
          auto newtree = new Tree();
          newtree->type = ATOMARRAY;
          if (fix->array_atom)
            newtree->array = &fix->array_atom[0][index1-1];
          newtree->nstride = fix->size_peratom_cols;
          treestack[ntreestack++] = newtree;
        } else print_var_error(FLERR,"Mismatched fix in variable formula",ivar);
      } else if (strncmp(word,"v_",2) == 0) {
        int ivar = find(word+2);
        if (ivar < 0)
          print_var_error(FLERR,fmt::format("Invalid variable reference {} in variable formula",word),
                          ivar);
        if (eval_in_progress[ivar])
          print_var_error(FLERR,"has a circular dependency",ivar);
        int nbracket;
        tagint index;
        if (str[i] != '[') nbracket = 0;
        else {
          nbracket = 1;
          ptr = &str[i];
          index = int_between_brackets(ptr,1);
          i = ptr-str+1;
        }
        if (nbracket == 0 && style[ivar] == INTERNAL) {
          value1 = dvalue[ivar];
          if (tree) {
            auto newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;
        } else if (nbracket == 0 && style[ivar] != ATOM &&
                   style[ivar] != VECTOR) {
          char *var = retrieve(word+2);
          if (var == nullptr)
            print_var_error(FLERR,"Invalid variable evaluation in variable formula",ivar);
          if (utils::is_double(var)) {
            if (tree) {
              auto newtree = new Tree();
              newtree->type = VALUE;
              newtree->value = atof(var);
              treestack[ntreestack++] = newtree;
            } else argstack[nargstack++] = atof(var);
          } else print_var_error(FLERR,"Non-numeric variable value in variable formula",ivar);
        } else if (nbracket == 0 && style[ivar] == ATOM) {
          if (tree == nullptr)
            print_var_error(FLERR,"Atom-style variable in equal-style variable formula",ivar);
          if (treetype == VECTOR)
            print_var_error(FLERR,"Atom-style variable in vector-style variable formula",ivar);
          Tree *newtree = nullptr;
          evaluate(data[ivar][0],&newtree,ivar);
          treestack[ntreestack++] = newtree;
        } else if (nbracket == 0 && style[ivar] == VECTOR) {
          if (tree == nullptr)
            print_var_error(FLERR,"Vector-style variable in equal-style variable formula",ivar);
          if (treetype == ATOM)
            print_var_error(FLERR,"Vector-style variable in atom-style variable formula",ivar);
          double *vec;
          int nvec = compute_vector(ivar,&vec);
          auto newtree = new Tree();
          newtree->type = VECTORARRAY;
          newtree->array = vec;
          newtree->nvector = nvec;
          newtree->nstride = 1;
          treestack[ntreestack++] = newtree;
        } else if (nbracket && style[ivar] == ATOM) {
          double *result;
          memory->create(result,atom->nlocal,"variable:result");
          compute_atom(ivar,0,result,1,0);
          peratom2global(1,nullptr,result,1,index,tree,treestack,ntreestack,argstack,nargstack);
          memory->destroy(result);
        } else if (nbracket && style[ivar] == VECTOR) {
          double *vec;
          int nvec = compute_vector(ivar,&vec);
          if (index <= 0 || index > nvec)
            print_var_error(FLERR,"Invalid index into vector-style variable",ivar);
          int m = index;
          if (tree) {
            auto newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = vec[m-1];
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = vec[m-1];
        } else print_var_error(FLERR,"Mismatched variable in variable formula",ivar);
      } else {
        if (str[i] == '(') {
          char *contents = nullptr;
          i = find_matching_paren(str,i,contents,ivar);
          i++;
          if (math_function(word,contents,tree,treestack,ntreestack,argstack,nargstack,ivar));
          else if (group_function(word,contents,tree,treestack,ntreestack,argstack,nargstack,ivar));
          else print_var_error(FLERR,fmt::format("Invalid math/group/special function '{}()' "
                                                 "in variable formula", word),ivar);
          delete[] contents;
        } else if (str[i] == '[') {
          if (domain->box_exist == 0)
            print_var_error(FLERR,"Variable evaluation before simulation box is defined",ivar);
          ptr = &str[i];
          tagint id = int_between_brackets(ptr,1);
          i = ptr-str+1;
          peratom2global(0,word,nullptr,0,id,tree,treestack,ntreestack,argstack,nargstack);
        } else if (is_atom_vector(word)) {
          if (domain->box_exist == 0)
            print_var_error(FLERR,"Variable evaluation before simulation box is defined",ivar);
          atom_vector(word,tree,treestack,ntreestack);
        } else if (constants.find(word) != constants.end()) {
          value1 = constants[word];
          if (tree) {
            auto newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;
        } else {
          if (domain->box_exist == 0)
            print_var_error(FLERR,"Variable evaluation before simulation box is defined",ivar);
          if (tree) {
            auto newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;
        }
      }
      delete[] word;
    } else if (strchr("+-*/^<>=!&|%\0",onechar)) {
      if (onechar == '+') op = ADD;
      else if (onechar == '-') op = SUBTRACT;
      else if (onechar == '*') op = MULTIPLY;
      else if (onechar == '/') op = DIVIDE;
      else if (onechar == '%') op = MODULO;
      else if (onechar == '^') op = CARAT;
      else if (onechar == '=') {
        if (str[i+1] != '=')
          print_var_error(FLERR,"Invalid syntax in variable formula",ivar);
        op = EQ;
        i++;
      } else if (onechar == '!') {
        if (str[i+1] == '=') {
          op = NE;
          i++;
        } else op = NOT;
      } else if (onechar == '<') {
        if (str[i+1] != '=') op = LT;
        else {
          op = LE;
          i++;
        }
      } else if (onechar == '>') {
        if (str[i+1] != '=') op = GT;
        else {
          op = GE;
          i++;
        }
      } else if (onechar == '&') {
        if (str[i+1] != '&')
          print_var_error(FLERR,"Invalid syntax in variable formula",ivar);
        op = AND;
        i++;
      } else if (onechar == '|') {
        if (str[i+1] == '|') op = OR;
        else if (str[i+1] == '^') op = XOR;
        else print_var_error(FLERR,"Invalid syntax in variable formula",ivar);
        i++;
      } else op = DONE;
      i++;
      if (op == SUBTRACT && expect == ARG) {
        opstack[nopstack++] = UNARY;
        continue;
      }
      if (op == NOT && expect == ARG) {
        opstack[nopstack++] = op;
        continue;
      }
      if (expect == ARG)
        print_var_error(FLERR,"Invalid syntax in variable formula",ivar);
      expect = ARG;
      while (nopstack && precedence[opstack[nopstack-1]] >= precedence[op]) {
        opprevious = opstack[--nopstack];
        if (tree) {
          auto newtree = new Tree();
          newtree->type = opprevious;
          if ((opprevious == UNARY) || (opprevious == NOT)) {
            newtree->first = treestack[--ntreestack];
          } else {
            newtree->second = treestack[--ntreestack];
            newtree->first = treestack[--ntreestack];
          }
          treestack[ntreestack++] = newtree;
        } else {
          value2 = argstack[--nargstack];
          if (opprevious != UNARY && opprevious != NOT)
            value1 = argstack[--nargstack];
          if (opprevious == ADD)
            argstack[nargstack++] = value1 + value2;
          else if (opprevious == SUBTRACT)
            argstack[nargstack++] = value1 - value2;
          else if (opprevious == MULTIPLY)
            argstack[nargstack++] = value1 * value2;
          else if (opprevious == DIVIDE) {
            if (value2 == 0.0)
              print_var_error(FLERR,"Divide by 0 in variable formula",ivar,0);
            argstack[nargstack++] = value1 / value2;
          } else if (opprevious == MODULO) {
            if (value2 == 0.0)
              print_var_error(FLERR,"Modulo 0 in variable formula",ivar,0);
            argstack[nargstack++] = fmod(value1,value2);
          } else if (opprevious == CARAT) {
            if (value2 == 0.0)
              argstack[nargstack++] = 1.0;
            else if ((value1 == 0.0) && (value2 < 0.0))
              print_var_error(FLERR,"Invalid power expression in variable formula",ivar,0);
            else argstack[nargstack++] = pow(value1,value2);
          } else if (opprevious == UNARY) {
            argstack[nargstack++] = -value2;
          } else if (opprevious == NOT) {
            if (value2 == 0.0) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          } else if (opprevious == EQ) {
            if (value1 == value2) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          } else if (opprevious == NE) {
            if (value1 != value2) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          } else if (opprevious == LT) {
            if (value1 < value2) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          } else if (opprevious == LE) {
            if (value1 <= value2) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          } else if (opprevious == GT) {
            if (value1 > value2) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          } else if (opprevious == GE) {
            if (value1 >= value2) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          } else if (opprevious == AND) {
            if (value1 != 0.0 && value2 != 0.0) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          } else if (opprevious == OR) {
            if (value1 != 0.0 || value2 != 0.0) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          } else if (opprevious == XOR) {
            if ((value1 == 0.0 && value2 != 0.0) ||
                (value1 != 0.0 && value2 == 0.0)) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          }
        }
      }
      if (op == DONE) break;
      opstack[nopstack++] = op;
    } else print_var_error(FLERR,"Invalid syntax in variable formula",ivar);
  }
  if (nopstack) print_var_error(FLERR,"Invalid syntax in variable formula",ivar);
  if (tree) {
    if (ntreestack != 1)
      print_var_error(FLERR,"Invalid syntax in variable formula",ivar);
    *tree = treestack[0];
    return 0.0;
  } else {
    if (nargstack != 1)
      print_var_error(FLERR,"Invalid syntax in variable formula",ivar);
    return argstack[0];
  }
}
double Variable::collapse_tree(Tree *tree)
{
  double arg1,arg2,arg3;
  if (tree->type == VALUE) return tree->value;
  if (tree->type == ATOMARRAY) return 0.0;
  if (tree->type == TYPEARRAY) return 0.0;
  if (tree->type == INTARRAY) return 0.0;
  if (tree->type == BIGINTARRAY) return 0.0;
  if (tree->type == VECTORARRAY) return 0.0;
  if (tree->type == ADD) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = arg1 + arg2;
    return tree->value;
  }
  if (tree->type == SUBTRACT) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = arg1 - arg2;
    return tree->value;
  }
  if (tree->type == MULTIPLY) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = arg1 * arg2;
    return tree->value;
  }
  if (tree->type == DIVIDE) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg2 == 0.0) error->one(FLERR,"Divide by 0 in variable formula");
    tree->value = arg1 / arg2;
    return tree->value;
  }
  if (tree->type == MODULO) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg2 == 0.0) error->one(FLERR,"Modulo 0 in variable formula");
    tree->value = fmod(arg1,arg2);
    return tree->value;
  }
  if (tree->type == CARAT) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg2 == 0.0) error->one(FLERR,"Power by 0 in variable formula");
    tree->value = pow(arg1,arg2);
    return tree->value;
  }
  if (tree->type == UNARY) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = -arg1;
    return tree->value;
  }
  if (tree->type == NOT) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 == 0.0) tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }
  if (tree->type == EQ) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 == arg2) tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }
  if (tree->type == NE) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 != arg2) tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }
  if (tree->type == LT) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 < arg2) tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }
  if (tree->type == LE) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 <= arg2) tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }
  if (tree->type == GT) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 > arg2) tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }
  if (tree->type == GE) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 >= arg2) tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }
  if (tree->type == AND) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 != 0.0 && arg2 != 0.0) tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }
  if (tree->type == OR) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 != 0.0 || arg2 != 0.0) tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }
  if (tree->type == XOR) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    if ((arg1 == 0.0 && arg2 != 0.0) || (arg1 != 0.0 && arg2 == 0.0))
      tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }
  if (tree->type == SQRT) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 < 0.0)
      error->one(FLERR,"Sqrt of negative value in variable formula");
    tree->value = sqrt(arg1);
    return tree->value;
  }
  if (tree->type == EXP) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = exp(arg1);
    return tree->value;
  }
  if (tree->type == LN) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 <= 0.0)
      error->one(FLERR,"Log of zero/negative value in variable formula");
    tree->value = log(arg1);
    return tree->value;
  }
  if (tree->type == LOG) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 <= 0.0)
      error->one(FLERR,"Log of zero/negative value in variable formula");
    tree->value = log10(arg1);
    return tree->value;
  }
  if (tree->type == ABS) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = fabs(arg1);
    return tree->value;
  }
  if (tree->type == SIN) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = sin(arg1);
    return tree->value;
  }
  if (tree->type == COS) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = cos(arg1);
    return tree->value;
  }
  if (tree->type == TAN) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = tan(arg1);
    return tree->value;
  }
  if (tree->type == ASIN) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 < -1.0 || arg1 > 1.0)
      error->one(FLERR,"Arcsin of invalid value in variable formula");
    tree->value = asin(arg1);
    return tree->value;
  }
  if (tree->type == ACOS) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 < -1.0 || arg1 > 1.0)
      error->one(FLERR,"Arccos of invalid value in variable formula");
    tree->value = acos(arg1);
    return tree->value;
  }
  if (tree->type == ATAN) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = atan(arg1);
    return tree->value;
  }
  if (tree->type == ATAN2) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = atan2(arg1,arg2);
    return tree->value;
  }
  if (tree->type == RANDOM) {
    collapse_tree(tree->first);
    collapse_tree(tree->second);
    if (randomatom == nullptr) {
      int seed = static_cast<int> (collapse_tree(tree->extra[0]));
      if (seed <= 0)
        error->one(FLERR,"Invalid math function in variable formula");
      randomatom = new RanMars(lmp,seed+me);
    }
    return 0.0;
  }
  if (tree->type == NORMAL) {
    collapse_tree(tree->first);
    double sigma = collapse_tree(tree->second);
    if (sigma < 0.0)
      error->one(FLERR,"Invalid math function in variable formula");
    if (randomatom == nullptr) {
      int seed = static_cast<int> (collapse_tree(tree->extra[0]));
      if (seed <= 0)
        error->one(FLERR,"Invalid math function in variable formula");
      randomatom = new RanMars(lmp,seed+me);
    }
    return 0.0;
  }
  if (tree->type == CEIL) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = ceil(arg1);
    return tree->value;
  }
  if (tree->type == FLOOR) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = floor(arg1);
    return tree->value;
  }
  if (tree->type == ROUND) {
    arg1 = collapse_tree(tree->first);
    if (tree->first->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = MYROUND(arg1);
    return tree->value;
  }
  if (tree->type == RAMP) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (update->whichflag == 0) {
      tree->value = arg1;
    } else {
      double delta = update->ntimestep - update->beginstep;
      if ((delta != 0.0) && (update->beginstep != update->endstep))
        delta /= update->endstep - update->beginstep;
      tree->value = arg1 + delta*(arg2-arg1);
    }
    return tree->value;
  }
  if (tree->type == STAGGER) {
    auto ivalue1 = static_cast<bigint> (collapse_tree(tree->first));
    auto ivalue2 = static_cast<bigint> (collapse_tree(tree->second));
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (ivalue1 <= 0 || ivalue2 <= 0 || ivalue1 <= ivalue2)
      error->one(FLERR,"Invalid math function in variable formula");
    bigint lower = update->ntimestep/ivalue1 * ivalue1;
    bigint delta = update->ntimestep - lower;
    if (delta < ivalue2) tree->value = lower+ivalue2;
    else tree->value = lower+ivalue1;
    return tree->value;
  }
  if (tree->type == LOGFREQ) {
    auto ivalue1 = static_cast<bigint> (collapse_tree(tree->first));
    auto ivalue2 = static_cast<bigint> (collapse_tree(tree->second));
    auto ivalue3 = static_cast<bigint> (collapse_tree(tree->extra[0]));
    if (tree->first->type != VALUE || tree->second->type != VALUE ||
        tree->extra[0]->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (ivalue1 <= 0 || ivalue2 <= 0 || ivalue3 <= 0 || ivalue2 >= ivalue3)
      error->one(FLERR,"Invalid math function in variable formula");
    if (update->ntimestep < ivalue1) tree->value = ivalue1;
    else {
      bigint lower = ivalue1;
      while (update->ntimestep >= ivalue3*lower) lower *= ivalue3;
      bigint multiple = update->ntimestep/lower;
      if (multiple < ivalue2) tree->value = (multiple+1)*lower;
      else tree->value = lower*ivalue3;
    }
    return tree->value;
  }
  if (tree->type == LOGFREQ2) {
    auto ivalue1 = static_cast<bigint> (collapse_tree(tree->first));
    auto ivalue2 = static_cast<bigint> (collapse_tree(tree->second));
    auto ivalue3 = static_cast<bigint> (collapse_tree(tree->extra[0]));
    if (tree->first->type != VALUE || tree->second->type != VALUE ||
        tree->extra[0]->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (ivalue1 <= 0 || ivalue2 <= 0 || ivalue3 <= 0 )
      error->all(FLERR,"Invalid math function in variable formula");
    if (update->ntimestep < ivalue1) tree->value = ivalue1;
    else {
      tree->value = ivalue1;
      double delta = ivalue1*(ivalue3-1.0)/ivalue2;
      bigint count = 0;
      while (update->ntimestep >= tree->value) {
        tree->value += delta;
        count++;
        if (count % ivalue2 == 0) delta *= ivalue3;
      }
    }
    tree->value = ceil(tree->value);
    return tree->value;
  }
  if (tree->type == LOGFREQ3) {
    auto ivalue1 = static_cast<bigint> (collapse_tree(tree->first));
    auto ivalue2 = static_cast<bigint> (collapse_tree(tree->second));
    auto ivalue3 = static_cast<bigint> (collapse_tree(tree->extra[0]));
    if (tree->first->type != VALUE || tree->second->type != VALUE ||
        tree->extra[0]->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (ivalue1 <= 0 || ivalue2 <= 1 || ivalue3 <= 0 ||
        ivalue3-ivalue1+1 < ivalue2 )
      error->all(FLERR,"Invalid math function in variable formula");
    if (update->ntimestep < ivalue1) tree->value = ivalue1;
    else {
      tree->value = ivalue1;
      double logsp = ivalue1;
      double factor = pow(((double)ivalue3)/ivalue1, 1.0/(ivalue2-1));
      bigint linsp = ivalue1;
      while (update->ntimestep >= (tree->value)) {
        logsp *= factor;
        linsp++;
        if (linsp > logsp) tree->value = linsp;
        else tree->value = ceil(logsp)-(((int)ceil(logsp)-1)/ivalue3);
      }
    }
    if (update->ntimestep > ivalue3)
      error->all(FLERR,"Calls to variable exceeded limit");
    return tree->value;
  }
  if (tree->type == STRIDE) {
    auto ivalue1 = static_cast<bigint> (collapse_tree(tree->first));
    auto ivalue2 = static_cast<bigint> (collapse_tree(tree->second));
    auto ivalue3 = static_cast<bigint> (collapse_tree(tree->extra[0]));
    if (tree->first->type != VALUE || tree->second->type != VALUE ||
        tree->extra[0]->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (ivalue1 < 0 || ivalue2 < 0 || ivalue3 <= 0 || ivalue1 > ivalue2)
      error->one(FLERR,"Invalid math function in variable formula");
    if (update->ntimestep < ivalue1) tree->value = ivalue1;
    else if (update->ntimestep < ivalue2) {
      bigint offset = update->ntimestep - ivalue1;
      tree->value = ivalue1 + (offset/ivalue3)*ivalue3 + ivalue3;
      if (tree->value > ivalue2) tree->value = (double) MAXBIGINT_DOUBLE;
    } else tree->value = (double) MAXBIGINT_DOUBLE;
    return tree->value;
  }
  if (tree->type == STRIDE2) {
    auto ivalue1 = static_cast<bigint> (collapse_tree(tree->first));
    auto ivalue2 = static_cast<bigint> (collapse_tree(tree->second));
    auto ivalue3 = static_cast<bigint> (collapse_tree(tree->extra[0]));
    auto ivalue4 = static_cast<bigint> (collapse_tree(tree->extra[1]));
    auto ivalue5 = static_cast<bigint> (collapse_tree(tree->extra[2]));
    auto ivalue6 = static_cast<bigint> (collapse_tree(tree->extra[3]));
    if (tree->first->type != VALUE || tree->second->type != VALUE ||
        tree->extra[0]->type != VALUE || tree->extra[1]->type != VALUE ||
        tree->extra[2]->type != VALUE || tree->extra[3]->type != VALUE)
      return 0.0;
    tree->type = VALUE;
    if (ivalue1 < 0 || ivalue2 < 0 || ivalue3 <= 0 || ivalue1 > ivalue2)
      error->one(FLERR,"Invalid math function in variable formula");
    if (ivalue4 < 0 || ivalue5 < 0 || ivalue6 <= 0 || ivalue4 > ivalue5)
      error->one(FLERR,"Invalid math function in variable formula");
    if (ivalue4 < ivalue1 || ivalue5 > ivalue2)
      error->one(FLERR,"Invalid math function in variable formula");
    bigint istep, offset;
    if (update->ntimestep < ivalue1) istep = ivalue1;
    else if (update->ntimestep < ivalue2) {
      if (update->ntimestep < ivalue4 || update->ntimestep > ivalue5) {
        offset = update->ntimestep - ivalue1;
        istep = ivalue1 + (offset/ivalue3)*ivalue3 + ivalue3;
        if (update->ntimestep < ivalue2 && istep > ivalue4)
          tree->value = ivalue4;
      } else {
        offset = update->ntimestep - ivalue4;
        istep = ivalue4 + (offset/ivalue6)*ivalue6 + ivalue6;
        if (istep > ivalue5) {
          offset = ivalue5 - ivalue1;
          istep = ivalue1 + (offset/ivalue3)*ivalue3 + ivalue3;
          if (istep > ivalue2) istep = MAXBIGINT_DOUBLE;
        }
      }
    } else istep = MAXBIGINT_DOUBLE;
    tree->value = (double)istep;
    return tree->value;
  }
  if (tree->type == VDISPLACE) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    if (tree->first->type != VALUE || tree->second->type != VALUE) return 0.0;
    tree->type = VALUE;
    double delta = update->ntimestep - update->beginstep;
    tree->value = arg1 + arg2*delta*update->dt;
    return tree->value;
  }
  if (tree->type == SWIGGLE) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    arg3 = collapse_tree(tree->extra[0]);
    if (tree->first->type != VALUE || tree->second->type != VALUE ||
        tree->extra[0]->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg3 == 0.0)
      error->one(FLERR,"Invalid swiggle(x,y,z) function in variable formula: z must be > 0");
    double delta = update->ntimestep - update->beginstep;
    double omega = 2.0*MY_PI/arg3;
    tree->value = arg1 + arg2*sin(omega*delta*update->dt);
    return tree->value;
  }
  if (tree->type == CWIGGLE) {
    arg1 = collapse_tree(tree->first);
    arg2 = collapse_tree(tree->second);
    arg3 = collapse_tree(tree->extra[0]);
    if (tree->first->type != VALUE || tree->second->type != VALUE ||
        tree->extra[0]->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg3 == 0.0)
      error->one(FLERR,"Invalid cwiggle(x,y,z) function in variable formula: z must be > 0");
    double delta = update->ntimestep - update->beginstep;
    double omega = 2.0*MY_PI/arg3;
    tree->value = arg1 + arg2*(1.0-cos(omega*delta*update->dt));
    return tree->value;
  }
  if (tree->type == GMASK) return 0.0;
  if (tree->type == RMASK) return 0.0;
  if (tree->type == GRMASK) return 0.0;
  return 0.0;
}
double Variable::eval_tree(Tree *tree, int i)
{
  double arg,arg1,arg2,arg3;
  if (tree->type == VALUE) return tree->value;
  if (tree->type == ATOMARRAY) return tree->array[i*tree->nstride];
  if (tree->type == TYPEARRAY) return tree->array[atom->type[i]];
  if (tree->type == INTARRAY) return (double) tree->iarray[i*tree->nstride];
  if (tree->type == BIGINTARRAY) return (double) tree->barray[i*tree->nstride];
  if (tree->type == VECTORARRAY) return tree->array[i*tree->nstride];
  if (tree->type == ADD)
    return eval_tree(tree->first,i) + eval_tree(tree->second,i);
  if (tree->type == SUBTRACT)
    return eval_tree(tree->first,i) - eval_tree(tree->second,i);
  if (tree->type == MULTIPLY)
    return eval_tree(tree->first,i) * eval_tree(tree->second,i);
  if (tree->type == DIVIDE) {
    double denom = eval_tree(tree->second,i);
    if (denom == 0.0) error->one(FLERR,"Divide by 0 in variable formula");
    return eval_tree(tree->first,i) / denom;
  }
  if (tree->type == MODULO) {
    double denom = eval_tree(tree->second,i);
    if (denom == 0.0) error->one(FLERR,"Modulo 0 in variable formula");
    return fmod(eval_tree(tree->first,i),denom);
  }
  if (tree->type == CARAT) {
    double exponent = eval_tree(tree->second,i);
    if (exponent == 0.0) error->one(FLERR,"Power by 0 in variable formula");
    return pow(eval_tree(tree->first,i),exponent);
  }
  if (tree->type == UNARY) return -eval_tree(tree->first,i);
  if (tree->type == NOT) {
    if (eval_tree(tree->first,i) == 0.0) return 1.0;
    else return 0.0;
  }
  if (tree->type == EQ) {
    if (eval_tree(tree->first,i) == eval_tree(tree->second,i)) return 1.0;
    else return 0.0;
  }
  if (tree->type == NE) {
    if (eval_tree(tree->first,i) != eval_tree(tree->second,i)) return 1.0;
    else return 0.0;
  }
  if (tree->type == LT) {
    if (eval_tree(tree->first,i) < eval_tree(tree->second,i)) return 1.0;
    else return 0.0;
  }
  if (tree->type == LE) {
    if (eval_tree(tree->first,i) <= eval_tree(tree->second,i)) return 1.0;
    else return 0.0;
  }
  if (tree->type == GT) {
    if (eval_tree(tree->first,i) > eval_tree(tree->second,i)) return 1.0;
    else return 0.0;
  }
  if (tree->type == GE) {
    if (eval_tree(tree->first,i) >= eval_tree(tree->second,i)) return 1.0;
    else return 0.0;
  }
  if (tree->type == AND) {
    if (eval_tree(tree->first,i) != 0.0 && eval_tree(tree->second,i) != 0.0)
      return 1.0;
    else return 0.0;
  }
  if (tree->type == OR) {
    if (eval_tree(tree->first,i) != 0.0 || eval_tree(tree->second,i) != 0.0)
      return 1.0;
    else return 0.0;
  }
  if (tree->type == XOR) {
    if ((eval_tree(tree->first,i) == 0.0 && eval_tree(tree->second,i) != 0.0)
        ||
        (eval_tree(tree->first,i) != 0.0 && eval_tree(tree->second,i) == 0.0))
      return 1.0;
    else return 0.0;
  }
  if (tree->type == SQRT) {
    arg1 = eval_tree(tree->first,i);
    if (arg1 < 0.0)
      error->one(FLERR,"Sqrt of negative value in variable formula");
    return sqrt(arg1);
  }
  if (tree->type == EXP)
    return exp(eval_tree(tree->first,i));
  if (tree->type == LN) {
    arg1 = eval_tree(tree->first,i);
    if (arg1 <= 0.0)
      error->one(FLERR,"Log of zero/negative value in variable formula");
    return log(arg1);
  }
  if (tree->type == LOG) {
    arg1 = eval_tree(tree->first,i);
    if (arg1 <= 0.0)
      error->one(FLERR,"Log of zero/negative value in variable formula");
    return log10(arg1);
  }
  if (tree->type == ABS)
    return fabs(eval_tree(tree->first,i));
  if (tree->type == SIN)
    return sin(eval_tree(tree->first,i));
  if (tree->type == COS)
    return cos(eval_tree(tree->first,i));
  if (tree->type == TAN)
    return tan(eval_tree(tree->first,i));
  if (tree->type == ASIN) {
    arg1 = eval_tree(tree->first,i);
    if (arg1 < -1.0 || arg1 > 1.0)
      error->one(FLERR,"Arcsin of invalid value in variable formula");
    return asin(arg1);
  }
  if (tree->type == ACOS) {
    arg1 = eval_tree(tree->first,i);
    if (arg1 < -1.0 || arg1 > 1.0)
      error->one(FLERR,"Arccos of invalid value in variable formula");
    return acos(arg1);
  }
  if (tree->type == ATAN)
    return atan(eval_tree(tree->first,i));
  if (tree->type == ATAN2)
    return atan2(eval_tree(tree->first,i),eval_tree(tree->second,i));
  if (tree->type == RANDOM) {
    double lower = eval_tree(tree->first,i);
    double upper = eval_tree(tree->second,i);
    if (randomatom == nullptr) {
      int seed = static_cast<int> (eval_tree(tree->extra[0],i));
      if (seed <= 0)
        error->one(FLERR,"Invalid math function in variable formula");
      randomatom = new RanMars(lmp,seed+me);
    }
    return randomatom->uniform()*(upper-lower)+lower;
  }
  if (tree->type == NORMAL) {
    double mu = eval_tree(tree->first,i);
    double sigma = eval_tree(tree->second,i);
    if (sigma < 0.0)
      error->one(FLERR,"Invalid math function in variable formula");
    if (randomatom == nullptr) {
      int seed = static_cast<int> (eval_tree(tree->extra[0],i));
      if (seed <= 0)
        error->one(FLERR,"Invalid math function in variable formula");
      randomatom = new RanMars(lmp,seed+me);
    }
    return mu + sigma*randomatom->gaussian();
  }
  if (tree->type == CEIL)
    return ceil(eval_tree(tree->first,i));
  if (tree->type == FLOOR)
    return floor(eval_tree(tree->first,i));
  if (tree->type == ROUND)
    return MYROUND(eval_tree(tree->first,i));
  if (tree->type == RAMP) {
    arg1 = eval_tree(tree->first,i);
    arg2 = eval_tree(tree->second,i);
    if (update->whichflag == 0) {
      arg = arg1;
    } else {
      double delta = update->ntimestep - update->beginstep;
      if ((delta != 0.0) && (update->beginstep != update->endstep))
        delta /= update->endstep - update->beginstep;
      arg = arg1 + delta*(arg2-arg1);
    }
    return arg;
  }
  if (tree->type == STAGGER) {
    auto ivalue1 = static_cast<bigint> (eval_tree(tree->first,i));
    auto ivalue2 = static_cast<bigint> (eval_tree(tree->second,i));
    if (ivalue1 <= 0 || ivalue2 <= 0 || ivalue1 <= ivalue2)
      error->one(FLERR,"Invalid math function in variable formula");
    bigint lower = update->ntimestep/ivalue1 * ivalue1;
    bigint delta = update->ntimestep - lower;
    if (delta < ivalue2) arg = lower+ivalue2;
    else arg = lower+ivalue1;
    return arg;
  }
  if (tree->type == LOGFREQ) {
    auto ivalue1 = static_cast<bigint> (eval_tree(tree->first,i));
    auto ivalue2 = static_cast<bigint> (eval_tree(tree->second,i));
    auto ivalue3 = static_cast<bigint> (eval_tree(tree->extra[0],i));
    if (ivalue1 <= 0 || ivalue2 <= 0 || ivalue3 <= 0 || ivalue2 >= ivalue3)
      error->one(FLERR,"Invalid math function in variable formula");
    if (update->ntimestep < ivalue1) arg = ivalue1;
    else {
      bigint lower = ivalue1;
      while (update->ntimestep >= ivalue3*lower) lower *= ivalue3;
      bigint multiple = update->ntimestep/lower;
      if (multiple < ivalue2) arg = (multiple+1)*lower;
      else arg = lower*ivalue3;
    }
    return arg;
  }
  if (tree->type == LOGFREQ2) {
    auto ivalue1 = static_cast<bigint> (eval_tree(tree->first,i));
    auto ivalue2 = static_cast<bigint> (eval_tree(tree->second,i));
    auto ivalue3 = static_cast<bigint> (eval_tree(tree->extra[0],i));
    if (ivalue1 <= 0 || ivalue2 <= 0 || ivalue3 <= 0 )
      error->all(FLERR,"Invalid math function in variable formula");
    if (update->ntimestep < ivalue1) arg = ivalue1;
    else {
      arg = ivalue1;
      double delta = ivalue1*(ivalue3-1.0)/ivalue2;
      bigint count = 0;
      while (update->ntimestep >= arg) {
        arg += delta;
        count++;
        if (count % ivalue2 == 0) delta *= ivalue3;
      }
    }
    arg = ceil(arg);
    return arg;
  }
  if (tree->type == STRIDE) {
    auto ivalue1 = static_cast<bigint> (eval_tree(tree->first,i));
    auto ivalue2 = static_cast<bigint> (eval_tree(tree->second,i));
    auto ivalue3 = static_cast<bigint> (eval_tree(tree->extra[0],i));
    if (ivalue1 < 0 || ivalue2 < 0 || ivalue3 <= 0 || ivalue1 > ivalue2)
      error->one(FLERR,"Invalid math function in variable formula");
    if (update->ntimestep < ivalue1) arg = ivalue1;
    else if (update->ntimestep < ivalue2) {
      bigint offset = update->ntimestep - ivalue1;
      arg = ivalue1 + (offset/ivalue3)*ivalue3 + ivalue3;
      if (arg > ivalue2) arg = (double) MAXBIGINT_DOUBLE;
    } else arg = (double) MAXBIGINT_DOUBLE;
    return arg;
  }
  if (tree->type == STRIDE2) {
    auto ivalue1 = static_cast<bigint> (eval_tree(tree->first,i));
    auto ivalue2 = static_cast<bigint> (eval_tree(tree->second,i));
    auto ivalue3 = static_cast<bigint> (eval_tree(tree->extra[0],i));
    auto ivalue4 = static_cast<bigint> (eval_tree(tree->extra[1],i));
    auto ivalue5 = static_cast<bigint> (eval_tree(tree->extra[2],i));
    auto ivalue6 = static_cast<bigint> (eval_tree(tree->extra[3],i));
    if (ivalue1 < 0 || ivalue2 < 0 || ivalue3 <= 0 || ivalue1 > ivalue2)
      error->one(FLERR,"Invalid math function in variable formula");
    if (ivalue4 < 0 || ivalue5 < 0 || ivalue6 <= 0 || ivalue4 > ivalue5)
      error->one(FLERR,"Invalid math function in variable formula");
    if (ivalue4 < ivalue1 || ivalue5 > ivalue2)
      error->one(FLERR,"Invalid math function in variable formula");
    bigint istep, offset;
    if (update->ntimestep < ivalue1) istep = ivalue1;
    else if (update->ntimestep < ivalue2) {
      if (update->ntimestep < ivalue4 || update->ntimestep > ivalue5) {
        offset = update->ntimestep - ivalue1;
        istep = ivalue1 + (offset/ivalue3)*ivalue3 + ivalue3;
        if (update->ntimestep < ivalue2 && istep > ivalue4)
          tree->value = ivalue4;
      } else {
        offset = update->ntimestep - ivalue4;
        istep = ivalue4 + (offset/ivalue6)*ivalue6 + ivalue6;
        if (istep > ivalue5) {
          offset = ivalue5 - ivalue1;
          istep = ivalue1 + (offset/ivalue3)*ivalue3 + ivalue3;
          if (istep > ivalue2) istep = MAXBIGINT_DOUBLE;
        }
      }
    } else istep = MAXBIGINT_DOUBLE;
    arg = istep;
    return arg;
  }
  if (tree->type == VDISPLACE) {
    arg1 = eval_tree(tree->first,i);
    arg2 = eval_tree(tree->second,i);
    double delta = update->ntimestep - update->beginstep;
    arg = arg1 + arg2*delta*update->dt;
    return arg;
  }
  if (tree->type == SWIGGLE) {
    arg1 = eval_tree(tree->first,i);
    arg2 = eval_tree(tree->second,i);
    arg3 = eval_tree(tree->extra[0],i);
    if (arg3 == 0.0)
      error->one(FLERR,"Invalid swiggle(x,y,z) function in variable formula: z must be > 0");
    double delta = update->ntimestep - update->beginstep;
    double omega = 2.0*MY_PI/arg3;
    arg = arg1 + arg2*sin(omega*delta*update->dt);
    return arg;
  }
  if (tree->type == CWIGGLE) {
    arg1 = eval_tree(tree->first,i);
    arg2 = eval_tree(tree->second,i);
    arg3 = eval_tree(tree->extra[0],i);
    if (arg3 == 0.0)
      error->one(FLERR,"Invalid cwiggle(x,y,z) function in variable formula: z must be > 0");
    double delta = update->ntimestep - update->beginstep;
    double omega = 2.0*MY_PI/arg3;
    arg = arg1 + arg2*(1.0-cos(omega*delta*update->dt));
    return arg;
  }
  if (tree->type == GMASK) {
    if (atom->mask[i] & tree->ivalue) return 1.0;
    else return 0.0;
  }
  if (tree->type == RMASK) {
    if (tree->region->match(atom->x[i][0], atom->x[i][1], atom->x[i][2])) return 1.0;
    else return 0.0;
  }
  if (tree->type == GRMASK) {
    if ((atom->mask[i] & tree->ivalue) &&
        (tree->region->match(atom->x[i][0], atom->x[i][1], atom->x[i][2]))) return 1.0;
    else return 0.0;
  }
  return 0.0;
}
int Variable::size_tree_vector(Tree *tree)
{
  int nsize = 0;
  if (tree->type == VECTORARRAY) nsize = tree->nvector;
  if (tree->first) nsize = compare_tree_vector(nsize, size_tree_vector(tree->first));
  if (tree->second) nsize = compare_tree_vector(nsize, size_tree_vector(tree->second));
  if (tree->nextra) {
    for (int i = 0; i < tree->nextra; i++)
      nsize = compare_tree_vector(nsize,size_tree_vector(tree->extra[i]));
  }
  return nsize;
}
int Variable::compare_tree_vector(int i, int j)
{
  if (i < 0 || j < 0) return -1;
  if (i == 0 || j == 0) return MAX(i,j);
  if (i != j) return -1;
  return i;
}
void Variable::free_tree(Tree *tree)
{
  if (tree->first) free_tree(tree->first);
  if (tree->second) free_tree(tree->second);
  if (tree->nextra) {
    for (int i = 0; i < tree->nextra; i++) free_tree(tree->extra[i]);
    delete[] tree->extra;
  }
  if (tree->selfalloc) memory->destroy(tree->array);
  delete tree;
}
int Variable::find_matching_paren(char *str, int i, char *&contents, int ivar)
{
  int istart = i;
  int ilevel = 0;
  while (true) {
    i++;
    if (!str[i]) break;
    if (str[i] == '(') ilevel++;
    else if (str[i] == ')' && ilevel) ilevel--;
    else if (str[i] == ')') break;
  }
  if (!str[i]) print_var_error(FLERR,"Invalid syntax in variable formula",ivar);
  int istop = i;
  int n = istop - istart - 1;
  delete[] contents;
  contents = new char[n+1];
  strncpy(contents,&str[istart+1],n);
  contents[n] = '\0';
  return istop;
}
tagint Variable::int_between_brackets(char *&ptr, int varallow)
{
  int varflag;
  tagint index;
  char *start = ++ptr;
  if (varallow && utils::strmatch(ptr,"^v_")) {
    varflag = 1;
    while (*ptr && *ptr != ']') {
      if (!isalnum(*ptr) && *ptr != '_')
        error->all(FLERR,"Variable name between brackets must be letters, numbers, or underscores");
      ptr++;
    }
  } else {
    varflag = 0;
    while (*ptr && *ptr != ']') {
      if (!isdigit(*ptr))
        error->all(FLERR,"Non digit character between brackets in variable");
      ptr++;
    }
  }
  if (*ptr != ']') error->all(FLERR,"Mismatched brackets in variable");
  if (ptr == start) error->all(FLERR,"Empty brackets in variable");
  *ptr = '\0';
  if (varflag) {
    char *id = start+2;
    int ivar = find(id);
    if (ivar < 0)
      error->all(FLERR,"Invalid variable name in variable formula");
    char *var = retrieve(id);
    if (var == nullptr)
      error->all(FLERR,"Invalid variable evaluation in variable formula");
    index = static_cast<tagint> (atof(var));
  } else index = ATOTAGINT(start);
  *ptr = ']';
  if (index == 0)
    error->all(FLERR,"Index between variable brackets must be positive");
  return index;
}
int Variable::math_function(char *word, char *contents, Tree **tree, Tree **treestack,
                            int &ntreestack, double *argstack, int &nargstack, int ivar)
{
  if (strcmp(word,"sqrt") != 0 && strcmp(word,"exp") &&
      strcmp(word,"ln") != 0 && strcmp(word,"log") != 0 &&
      strcmp(word,"abs") != 0 &&
      strcmp(word,"sin") != 0 && strcmp(word,"cos") != 0 &&
      strcmp(word,"tan") != 0 && strcmp(word,"asin") != 0 &&
      strcmp(word,"acos") != 0 && strcmp(word,"atan") != 0 &&
      strcmp(word,"atan2") != 0 && strcmp(word,"random") != 0 &&
      strcmp(word,"normal") != 0 && strcmp(word,"ceil") != 0 &&
      strcmp(word,"floor") != 0 && strcmp(word,"round") != 0 &&
      strcmp(word,"ramp") != 0 && strcmp(word,"stagger") != 0 &&
      strcmp(word,"logfreq") != 0 && strcmp(word,"logfreq2") != 0 &&
      strcmp(word,"logfreq3") != 0 && strcmp(word,"stride") != 0 &&
      strcmp(word,"stride2") != 0 && strcmp(word,"vdisplace") != 0 &&
      strcmp(word,"swiggle") != 0 && strcmp(word,"cwiggle") != 0)
    return 0;
  char *args[MAXFUNCARG];
  int narg = parse_args(contents,args);
  Tree *newtree = nullptr;
  double value1,value2;
  double values[MAXFUNCARG-2];
  if (tree) {
    newtree = new Tree();
    Tree *argtree = nullptr;
    evaluate(args[0],&argtree,ivar);
    newtree->first = argtree;
    if (narg > 1) {
      evaluate(args[1],&argtree,ivar);
      newtree->second = argtree;
      if (narg > 2) {
        newtree->nextra = narg-2;
        newtree->extra = new Tree*[narg-2];
        for (int i = 2; i < narg; i++) {
          evaluate(args[i],&argtree,ivar);
          newtree->extra[i-2] = argtree;
        }
      }
    }
    treestack[ntreestack++] = newtree;
  } else {
    value1 = evaluate(args[0],nullptr,ivar);
    if (narg > 1) {
      value2 = evaluate(args[1],nullptr,ivar);
      if (narg > 2) {
        for (int i = 2; i < narg; i++)
          values[i-2] = evaluate(args[i],nullptr,ivar);
      }
    }
  }
  if (strcmp(word,"sqrt") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = SQRT;
    else {
      if (value1 < 0.0)
        print_var_error(FLERR,"Sqrt of negative value in variable formula",ivar,0);
      argstack[nargstack++] = sqrt(value1);
    }
  } else if (strcmp(word,"exp") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = EXP;
    else argstack[nargstack++] = exp(value1);
  } else if (strcmp(word,"ln") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = LN;
    else {
      if (value1 <= 0.0)
        print_var_error(FLERR,"Log of zero/negative value in variable formula",ivar,0);
      argstack[nargstack++] = log(value1);
    }
  } else if (strcmp(word,"log") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = LOG;
    else {
      if (value1 <= 0.0)
        print_var_error(FLERR,"Log of zero/negative value in variable formula",ivar,0);
      argstack[nargstack++] = log10(value1);
    }
  } else if (strcmp(word,"abs") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = ABS;
    else argstack[nargstack++] = fabs(value1);
  } else if (strcmp(word,"sin") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = SIN;
    else argstack[nargstack++] = sin(value1);
  } else if (strcmp(word,"cos") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = COS;
    else argstack[nargstack++] = cos(value1);
  } else if (strcmp(word,"tan") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = TAN;
    else argstack[nargstack++] = tan(value1);
  } else if (strcmp(word,"asin") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = ASIN;
    else {
      if (value1 < -1.0 || value1 > 1.0)
        print_var_error(FLERR,"Arcsin of invalid value in variable formula",ivar,0);
      argstack[nargstack++] = asin(value1);
    }
  } else if (strcmp(word,"acos") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = ACOS;
    else {
      if (value1 < -1.0 || value1 > 1.0)
        print_var_error(FLERR,"Arccos of invalid value in variable formula",ivar,0);
      argstack[nargstack++] = acos(value1);
    }
  } else if (strcmp(word,"atan") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = ATAN;
    else argstack[nargstack++] = atan(value1);
  } else if (strcmp(word,"atan2") == 0) {
    if (narg != 2)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = ATAN2;
    else argstack[nargstack++] = atan2(value1,value2);
  } else if (strcmp(word,"random") == 0) {
    if (narg != 3)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = RANDOM;
    else {
      if (randomequal == nullptr) {
        int seed = static_cast<int> (values[0]);
        if (seed <= 0)
          print_var_error(FLERR,"Invalid math function in variable formula",ivar);
        randomequal = new RanMars(lmp,seed);
      }
      argstack[nargstack++] = randomequal->uniform()*(value2-value1) + value1;
    }
  } else if (strcmp(word,"normal") == 0) {
    if (narg != 3)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = NORMAL;
    else {
      if (value2 < 0.0)
        print_var_error(FLERR,"Invalid math function in variable formula",ivar);
      if (randomequal == nullptr) {
        int seed = static_cast<int> (values[0]);
        if (seed <= 0)
          print_var_error(FLERR,"Invalid math function in variable formula",ivar);
        randomequal = new RanMars(lmp,seed);
      }
      argstack[nargstack++] = value1 + value2*randomequal->gaussian();
    }
  } else if (strcmp(word,"ceil") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = CEIL;
    else argstack[nargstack++] = ceil(value1);
  } else if (strcmp(word,"floor") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = FLOOR;
    else argstack[nargstack++] = floor(value1);
  } else if (strcmp(word,"round") == 0) {
    if (narg != 1)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = ROUND;
    else argstack[nargstack++] = MYROUND(value1);
  } else if (strcmp(word,"ramp") == 0) {
    if (narg != 2)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = RAMP;
    else {
      if (update->whichflag == 0) {
        argstack[nargstack++] = value1;
      } else {
        double delta = update->ntimestep - update->beginstep;
        if ((delta != 0.0) && (update->beginstep != update->endstep))
          delta /= update->endstep - update->beginstep;
        double value = value1 + delta*(value2-value1);
        argstack[nargstack++] = value;
      }
    }
  } else if (strcmp(word,"stagger") == 0) {
    if (narg != 2)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = STAGGER;
    else {
      auto ivalue1 = static_cast<bigint> (value1);
      auto ivalue2 = static_cast<bigint> (value2);
      if (ivalue1 <= 0 || ivalue2 <= 0 || ivalue1 <= ivalue2)
        print_var_error(FLERR,"Invalid math function in variable formula",ivar);
      bigint lower = update->ntimestep/ivalue1 * ivalue1;
      bigint delta = update->ntimestep - lower;
      double value;
      if (delta < ivalue2) value = lower+ivalue2;
      else value = lower+ivalue1;
      argstack[nargstack++] = value;
    }
  } else if (strcmp(word,"logfreq") == 0) {
    if (narg != 3)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = LOGFREQ;
    else {
      auto ivalue1 = static_cast<bigint> (value1);
      auto ivalue2 = static_cast<bigint> (value2);
      auto ivalue3 = static_cast<bigint> (values[0]);
      if (ivalue1 <= 0 || ivalue2 <= 0 || ivalue3 <= 0 || ivalue2 >= ivalue3)
        print_var_error(FLERR,"Invalid math function in variable formula",ivar);
      double value;
      if (update->ntimestep < ivalue1) value = ivalue1;
      else {
        bigint lower = ivalue1;
        while (update->ntimestep >= ivalue3*lower) lower *= ivalue3;
        bigint multiple = update->ntimestep/lower;
        if (multiple < ivalue2) value = (multiple+1)*lower;
        else value = lower*ivalue3;
      }
      argstack[nargstack++] = value;
    }
  } else if (strcmp(word,"logfreq2") == 0) {
    if (narg != 3)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = LOGFREQ2;
    else {
      auto ivalue1 = static_cast<bigint> (value1);
      auto ivalue2 = static_cast<bigint> (value2);
      auto ivalue3 = static_cast<bigint> (values[0]);
      if (ivalue1 <= 0 || ivalue2 <= 0 || ivalue3 <= 0 )
        print_var_error(FLERR,"Invalid math function in variable formula",ivar);
      double value;
      if (update->ntimestep < ivalue1) value = ivalue1;
      else {
        value = ivalue1;
        double delta = ivalue1*(ivalue3-1.0)/ivalue2;
        bigint count = 0;
        while (update->ntimestep >= value) {
          value += delta;
          count++;
          if (count % ivalue2 == 0) delta *= ivalue3;
        }
      }
      argstack[nargstack++] = ceil(value);
    }
  } else if (strcmp(word,"logfreq3") == 0) {
    if (narg != 3)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = LOGFREQ3;
    else {
      auto ivalue1 = static_cast<bigint> (value1);
      auto ivalue2 = static_cast<bigint> (value2);
      auto ivalue3 = static_cast<bigint> (values[0]);
      if (ivalue1 <= 0 || ivalue2 <= 1 || ivalue3 <= 0 ||
          ivalue3-ivalue1+1 < ivalue2 )
        print_var_error(FLERR,"Invalid math function in variable formula",ivar);
      double value;
      if (update->ntimestep < ivalue1) value = ivalue1;
      else {
        value = ivalue1;
        double logsp = ivalue1;
        double factor = pow(((double)ivalue3)/ivalue1, 1.0/(ivalue2-1));
        bigint linsp = ivalue1;
        while (update->ntimestep >= value) {
          logsp *= factor;
          linsp++;
          if (linsp > logsp) value = linsp;
          else value = ceil(logsp)-(((bigint)ceil(logsp)-1)/ivalue3);
        }
      }
      if (update->ntimestep > ivalue3)
        error->all(FLERR,"Calls to variable exceeded limit");
      argstack[nargstack++] = value;
    }
  } else if (strcmp(word,"stride") == 0) {
    if (narg != 3)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = STRIDE;
    else {
      auto ivalue1 = static_cast<bigint> (value1);
      auto ivalue2 = static_cast<bigint> (value2);
      auto ivalue3 = static_cast<bigint> (values[0]);
      if (ivalue1 < 0 || ivalue2 < 0 || ivalue3 <= 0 || ivalue1 > ivalue2)
        error->one(FLERR,"Invalid math function in variable formula");
      double value;
      if (update->ntimestep < ivalue1) value = ivalue1;
      else if (update->ntimestep < ivalue2) {
        bigint offset = update->ntimestep - ivalue1;
        value = ivalue1 + (offset/ivalue3)*ivalue3 + ivalue3;
        if (value > ivalue2) value = (double) MAXBIGINT_DOUBLE;
      } else value = (double) MAXBIGINT_DOUBLE;
      argstack[nargstack++] = value;
    }
  } else if (strcmp(word,"stride2") == 0) {
    if (narg != 6)
      print_var_error(FLERR,"Invalid math function in variable formula",ivar);
    if (tree) newtree->type = STRIDE2;
    else {
      auto ivalue1 = static_cast<bigint> (value1);
      auto ivalue2 = static_cast<bigint> (value2);
      auto ivalue3 = static_cast<bigint> (values[0]);
      auto ivalue4 = static_cast<bigint> (values[1]);
      auto ivalue5 = static_cast<bigint> (values[2]);
      auto ivalue6 = static_cast<bigint> (values[3]);
      if (ivalue1 < 0 || ivalue2 < 0 || ivalue3 <= 0 || ivalue1 > ivalue2)
        error->one(FLERR,"Invalid math function in variable formula");
      if (ivalue4 < 0 || ivalue5 < 0 || ivalue6 <= 0 || ivalue4 > ivalue5)
        error->one(FLERR,"Invalid math function in variable formula");
      if (ivalue4 < ivalue1 || ivalue5 > ivalue2)
        error->one(FLERR,"Invalid math function in variable formula");
      bigint istep, offset;
      if (update->ntimestep < ivalue1) istep = ivalue1;
      else if (update->ntimestep < ivalue2) {
        if (update->ntimestep < ivalue4 || update->ntimestep > ivalue5) {
          offset = update->ntimestep - ivalue1;
          istep = ivalue1 + (offset/ivalue3)*ivalue3 + ivalue3;
          if (update->ntimestep < ivalue4 && istep > ivalue4) istep = ivalue4;
        } else {
          offset = update->ntimestep - ivalue4;
          istep = ivalue4 + (offset/ivalue6)*ivalue6 + ivalue6;
          if (istep > ivalue5) {
            offset = ivalue5 - ivalue1;
            istep = ivalue1 + (offset/ivalue3)*ivalue3 + ivalue3;
            if (istep > ivalue2) istep = MAXBIGINT_DOUBLE;
          }
        }
      } else istep = MAXBIGINT_DOUBLE;
      double value = istep;
      argstack[nargstack++] = value;
    }
  } else if (strcmp(word,"vdisplace") == 0) {
    if (narg != 2)
      print_var_error(FLERR,"Invalid vdisplace function in variable formula: must have 2 arguments",ivar);
    if (modify->get_fix_by_style("dt/reset").size() > 0)
      print_var_error(FLERR,"Must not use vdisplace(x,y) function with fix dt/reset",ivar);
    if (tree) newtree->type = VDISPLACE;
    else {
      double delta = update->ntimestep - update->beginstep;
      double value = value1 + value2*delta*update->dt;
      argstack[nargstack++] = value;
    }
  } else if (strcmp(word,"swiggle") == 0) {
    if (narg != 3)
      print_var_error(FLERR,"Invalid swiggle function in variable formula: must have 3 arguments",ivar);
    if (modify->get_fix_by_style("dt/reset").size() > 0)
      print_var_error(FLERR,"Must not use swiggle(x,y,z) function with fix dt/reset",ivar);
    if (tree) newtree->type = SWIGGLE;
    else {
      if (values[0] == 0.0)
        print_var_error(FLERR,"Invalid swiggle(x,y,z) function in variable formula: z must be > 0",ivar);
      double delta = update->ntimestep - update->beginstep;
      double omega = 2.0*MY_PI/values[0];
      double value = value1 + value2*sin(omega*delta*update->dt);
      argstack[nargstack++] = value;
    }
  } else if (strcmp(word,"cwiggle") == 0) {
    if (narg != 3)
      print_var_error(FLERR,"Invalid cwiggle function in variable formula: must have 3 arguments",ivar);
    if (modify->get_fix_by_style("dt/reset").size() > 0)
      print_var_error(FLERR,"Must not use cwiggle(x,y,z) function with fix dt/reset",ivar);
    if (tree) newtree->type = CWIGGLE;
    else {
      if (values[0] == 0.0)
        print_var_error(FLERR,"Invalid cwiggle(x,y,z) function in variable formula: z must be > 0",ivar);
      double delta = update->ntimestep - update->beginstep;
      double omega = 2.0*MY_PI/values[0];
      double value = value1 + value2*(1.0-cos(omega*delta*update->dt));
      argstack[nargstack++] = value;
    }
  }
  for (int i = 0; i < narg; i++) delete[] args[i];
  return 1;
}
int Variable::group_function(char *word, char *contents, Tree **tree, Tree **treestack,
                             int &ntreestack, double *argstack, int &nargstack, int ivar)
{
  if (strcmp(word,"count") != 0 && strcmp(word,"mass") &&
      strcmp(word,"charge") != 0 && strcmp(word,"xcm") != 0 &&
      strcmp(word,"vcm") != 0 && strcmp(word,"fcm") != 0 &&
      strcmp(word,"bound") != 0 && strcmp(word,"gyration") != 0 &&
      strcmp(word,"ke") != 0 && strcmp(word,"angmom") != 0 &&
      strcmp(word,"torque") != 0 && strcmp(word,"inertia") != 0 &&
      strcmp(word,"omega") != 0)
    return 0;
  char *args[MAXFUNCARG];
  int narg = parse_args(contents,args);
  int igroup = group->find(args[0]);
  if (igroup == -1) {
    const auto errmesg = fmt::format("Group {} in variable formula does not exist", args[0]);
    print_var_error(FLERR, errmesg, ivar);
  }
  double value = 0.0;
  const auto group_errmesg = fmt::format("Invalid {}() function in variable formula", word);
  if (strcmp(word,"count") == 0) {
    if (narg == 1) value = group->count(igroup);
    else if (narg == 2)
      value = group->count(igroup,region_function(args[1],ivar));
    else print_var_error(FLERR,group_errmesg,ivar);
  } else if (strcmp(word,"mass") == 0) {
    if (narg == 1) value = group->mass(igroup);
    else if (narg == 2) value = group->mass(igroup,region_function(args[1],ivar));
    else print_var_error(FLERR,group_errmesg,ivar);
  } else if (strcmp(word,"charge") == 0) {
    if (narg == 1) value = group->charge(igroup);
    else if (narg == 2)
      value = group->charge(igroup,region_function(args[1],ivar));
    else print_var_error(FLERR,group_errmesg,ivar);
  } else if (strcmp(word,"xcm") == 0) {
    atom->check_mass(FLERR);
    double xcm[3];
    if (narg == 2) {
      double masstotal = group->mass(igroup);
      group->xcm(igroup,masstotal,xcm);
    } else if (narg == 3) {
      auto region = region_function(args[2],ivar);
      double masstotal = group->mass(igroup,region);
      group->xcm(igroup,masstotal,xcm,region);
    } else print_var_error(FLERR,group_errmesg,ivar);
    if (strcmp(args[1],"x") == 0) value = xcm[0];
    else if (strcmp(args[1],"y") == 0) value = xcm[1];
    else if (strcmp(args[1],"z") == 0) value = xcm[2];
    else print_var_error(FLERR,group_errmesg,ivar);
  } else if (strcmp(word,"vcm") == 0) {
    atom->check_mass(FLERR);
    double vcm[3];
    if (narg == 2) {
      double masstotal = group->mass(igroup);
      group->vcm(igroup,masstotal,vcm);
    } else if (narg == 3) {
      auto region = region_function(args[2],ivar);
      double masstotal = group->mass(igroup,region);
      group->vcm(igroup,masstotal,vcm,region);
    } else print_var_error(FLERR,group_errmesg,ivar);
    if (strcmp(args[1],"x") == 0) value = vcm[0];
    else if (strcmp(args[1],"y") == 0) value = vcm[1];
    else if (strcmp(args[1],"z") == 0) value = vcm[2];
    else print_var_error(FLERR,group_errmesg,ivar);
  } else if (strcmp(word,"fcm") == 0) {
    double fcm[3];
    if (narg == 2) group->fcm(igroup,fcm);
    else if (narg == 3) group->fcm(igroup,fcm,region_function(args[2],ivar));
    else print_var_error(FLERR,group_errmesg,ivar);
    if (strcmp(args[1],"x") == 0) value = fcm[0];
    else if (strcmp(args[1],"y") == 0) value = fcm[1];
    else if (strcmp(args[1],"z") == 0) value = fcm[2];
    else print_var_error(FLERR,group_errmesg,ivar);
  } else if (strcmp(word,"bound") == 0) {
    double minmax[6];
    if (narg == 2) group->bounds(igroup,minmax);
    else if (narg == 3)
      group->bounds(igroup,minmax,region_function(args[2],ivar));
    else print_var_error(FLERR,group_errmesg,ivar);
    if (strcmp(args[1],"xmin") == 0) value = minmax[0];
    else if (strcmp(args[1],"xmax") == 0) value = minmax[1];
    else if (strcmp(args[1],"ymin") == 0) value = minmax[2];
    else if (strcmp(args[1],"ymax") == 0) value = minmax[3];
    else if (strcmp(args[1],"zmin") == 0) value = minmax[4];
    else if (strcmp(args[1],"zmax") == 0) value = minmax[5];
    else print_var_error(FLERR,group_errmesg,ivar);
  } else if (strcmp(word,"gyration") == 0) {
    atom->check_mass(FLERR);
    double xcm[3];
    if (narg == 1) {
      double masstotal = group->mass(igroup);
      group->xcm(igroup,masstotal,xcm);
      value = group->gyration(igroup,masstotal,xcm);
    } else if (narg == 2) {
      auto region = region_function(args[1],ivar);
      double masstotal = group->mass(igroup,region);
      group->xcm(igroup,masstotal,xcm,region);
      value = group->gyration(igroup,masstotal,xcm,region);
    } else print_var_error(FLERR,group_errmesg,ivar);
  } else if (strcmp(word,"ke") == 0) {
    if (narg == 1) value = group->ke(igroup);
    else if (narg == 2) value = group->ke(igroup,region_function(args[1],ivar));
    else print_var_error(FLERR,group_errmesg,ivar);
  } else if (strcmp(word,"angmom") == 0) {
    atom->check_mass(FLERR);
    double xcm[3],lmom[3];
    if (narg == 2) {
      double masstotal = group->mass(igroup);
      group->xcm(igroup,masstotal,xcm);
      group->angmom(igroup,xcm,lmom);
    } else if (narg == 3) {
      auto region = region_function(args[2],ivar);
      double masstotal = group->mass(igroup,region);
      group->xcm(igroup,masstotal,xcm,region);
      group->angmom(igroup,xcm,lmom,region);
    } else print_var_error(FLERR,group_errmesg,ivar);
    if (strcmp(args[1],"x") == 0) value = lmom[0];
    else if (strcmp(args[1],"y") == 0) value = lmom[1];
    else if (strcmp(args[1],"z") == 0) value = lmom[2];
    else print_var_error(FLERR,group_errmesg,ivar);
  } else if (strcmp(word,"torque") == 0) {
    atom->check_mass(FLERR);
    double xcm[3],tq[3];
    if (narg == 2) {
      double masstotal = group->mass(igroup);
      group->xcm(igroup,masstotal,xcm);
      group->torque(igroup,xcm,tq);
    } else if (narg == 3) {
      auto region = region_function(args[2],ivar);
      double masstotal = group->mass(igroup,region);
      group->xcm(igroup,masstotal,xcm,region);
      group->torque(igroup,xcm,tq,region);
    } else print_var_error(FLERR,group_errmesg,ivar);
    if (strcmp(args[1],"x") == 0) value = tq[0];
    else if (strcmp(args[1],"y") == 0) value = tq[1];
    else if (strcmp(args[1],"z") == 0) value = tq[2];
    else print_var_error(FLERR,group_errmesg,ivar);
  } else if (strcmp(word,"inertia") == 0) {
    atom->check_mass(FLERR);
    double xcm[3],inertia[3][3];
    if (narg == 2) {
      double masstotal = group->mass(igroup);
      group->xcm(igroup,masstotal,xcm);
      group->inertia(igroup,xcm,inertia);
    } else if (narg == 3) {
      auto region = region_function(args[2],ivar);
      double masstotal = group->mass(igroup,region);
      group->xcm(igroup,masstotal,xcm,region);
      group->inertia(igroup,xcm,inertia,region);
    } else print_var_error(FLERR,group_errmesg,ivar);
    if (strcmp(args[1],"xx") == 0) value = inertia[0][0];
    else if (strcmp(args[1],"yy") == 0) value = inertia[1][1];
    else if (strcmp(args[1],"zz") == 0) value = inertia[2][2];
    else if (strcmp(args[1],"xy") == 0) value = inertia[0][1];
    else if (strcmp(args[1],"yz") == 0) value = inertia[1][2];
    else if (strcmp(args[1],"xz") == 0) value = inertia[0][2];
    else print_var_error(FLERR,group_errmesg,ivar);
  } else if (strcmp(word,"omega") == 0) {
    atom->check_mass(FLERR);
    double xcm[3],angmom[3],inertia[3][3],omega[3];
    if (narg == 2) {
      double masstotal = group->mass(igroup);
      group->xcm(igroup,masstotal,xcm);
      group->angmom(igroup,xcm,angmom);
      group->inertia(igroup,xcm,inertia);
      group->omega(angmom,inertia,omega);
    } else if (narg == 3) {
      auto region = region_function(args[2],ivar);
      double masstotal = group->mass(igroup,region);
      group->xcm(igroup,masstotal,xcm,region);
      group->angmom(igroup,xcm,angmom,region);
      group->inertia(igroup,xcm,inertia,region);
      group->omega(angmom,inertia,omega);
    } else print_var_error(FLERR,group_errmesg,ivar);
    if (strcmp(args[1],"x") == 0) value = omega[0];
    else if (strcmp(args[1],"y") == 0) value = omega[1];
    else if (strcmp(args[1],"z") == 0) value = omega[2];
    else print_var_error(FLERR,group_errmesg,ivar);
  }
  for (int i = 0; i < narg; i++) delete[] args[i];
  if (tree) {
    auto newtree = new Tree();
    newtree->type = VALUE;
    newtree->value = value;
    treestack[ntreestack++] = newtree;
  } else argstack[nargstack++] = value;
  return 1;
}
Region *Variable::region_function(char *id, int ivar)
{
  auto region = domain->get_region_by_id(id);
  if (!region)
    print_var_error(FLERR, fmt::format("Region {} in variable formula does not exist", id), ivar);
  region->init();
  return region;
}
void Variable::peratom2global(int flag, char *word, double *vector, int nstride, tagint id, Tree **tree,
                              Tree **treestack, int &ntreestack, double *argstack, int &nargstack)
{
  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR, "Indexed per-atom vector in variable formula without atom map");
  if (id > atom->map_tag_max)
    error->all(FLERR,"Variable atom ID is too large");
  int index = atom->map(id);
  double mine;
  if (index >= 0 && index < atom->nlocal) {
    if (flag == 0) {
      if (strcmp(word,"id") == 0) mine = atom->tag[index];
      else if (strcmp(word,"mass") == 0) {
        if (atom->rmass) mine = atom->rmass[index];
        else mine = atom->mass[atom->type[index]];
      }
      else if (strcmp(word,"type") == 0) mine = atom->type[index];
      else if (strcmp(word,"x") == 0) mine = atom->x[index][0];
      else if (strcmp(word,"y") == 0) mine = atom->x[index][1];
      else if (strcmp(word,"z") == 0) mine = atom->x[index][2];
      else if (strcmp(word,"vx") == 0) mine = atom->v[index][0];
      else if (strcmp(word,"vy") == 0) mine = atom->v[index][1];
      else if (strcmp(word,"vz") == 0) mine = atom->v[index][2];
      else if (strcmp(word,"fx") == 0) mine = atom->f[index][0];
      else if (strcmp(word,"fy") == 0) mine = atom->f[index][1];
      else if (strcmp(word,"fz") == 0) mine = atom->f[index][2];
      else if (strcmp(word,"q") == 0) {
        if (!atom->q_flag)
          error->one(FLERR,"Variable uses atom property that isn't allocated");
        mine = atom->q[index];
      }
      else error->one(FLERR,"Invalid atom vector {} in variable formula", word);
    } else mine = vector[index*nstride];
  } else mine = 0.0;
  double value;
  MPI_Allreduce(&mine,&value,1,MPI_DOUBLE,MPI_SUM,world);
  if (tree) {
    auto newtree = new Tree();
    newtree->type = VALUE;
    newtree->value = value;
    treestack[ntreestack++] = newtree;
  } else argstack[nargstack++] = value;
}
int Variable::is_atom_vector(char *word)
{
  if (strcmp(word,"id") == 0) return 1;
  if (strcmp(word,"mass") == 0) return 1;
  if (strcmp(word,"type") == 0) return 1;
  if (strcmp(word,"mol") == 0) return 1;
  if (strcmp(word,"x") == 0) return 1;
  if (strcmp(word,"y") == 0) return 1;
  if (strcmp(word,"z") == 0) return 1;
  if (strcmp(word,"vx") == 0) return 1;
  if (strcmp(word,"vy") == 0) return 1;
  if (strcmp(word,"vz") == 0) return 1;
  if (strcmp(word,"fx") == 0) return 1;
  if (strcmp(word,"fy") == 0) return 1;
  if (strcmp(word,"fz") == 0) return 1;
  if (strcmp(word,"q") == 0) return 1;
  return 0;
}
void Variable::atom_vector(char *word, Tree **tree, Tree **treestack, int &ntreestack)
{
  if (tree == nullptr)
    error->all(FLERR,"Atom vector in equal-style variable formula");
  auto newtree = new Tree();
  newtree->type = ATOMARRAY;
  newtree->nstride = 3;
  treestack[ntreestack++] = newtree;
  if (strcmp(word,"id") == 0) {
    if (sizeof(tagint) == sizeof(smallint)) {
      newtree->type = INTARRAY;
      newtree->iarray = (int *) atom->tag;
    } else {
      newtree->type = BIGINTARRAY;
      newtree->barray = (bigint *) atom->tag;
    }
    newtree->nstride = 1;
  } else if (strcmp(word,"mass") == 0) {
    if (atom->rmass) {
      newtree->nstride = 1;
      newtree->array = atom->rmass;
    } else {
      newtree->type = TYPEARRAY;
      newtree->array = atom->mass;
    }
  } else if (strcmp(word,"type") == 0) {
    newtree->type = INTARRAY;
    newtree->nstride = 1;
    newtree->iarray = atom->type;
  }
  else if (strcmp(word,"x") == 0) newtree->array = &atom->x[0][0];
  else if (strcmp(word,"y") == 0) newtree->array = &atom->x[0][1];
  else if (strcmp(word,"z") == 0) newtree->array = &atom->x[0][2];
  else if (strcmp(word,"vx") == 0) newtree->array = &atom->v[0][0];
  else if (strcmp(word,"vy") == 0) newtree->array = &atom->v[0][1];
  else if (strcmp(word,"vz") == 0) newtree->array = &atom->v[0][2];
  else if (strcmp(word,"fx") == 0) newtree->array = &atom->f[0][0];
  else if (strcmp(word,"fy") == 0) newtree->array = &atom->f[0][1];
  else if (strcmp(word,"fz") == 0) newtree->array = &atom->f[0][2];
  else if (strcmp(word,"q") == 0) {
    newtree->nstride = 1;
    newtree->array = atom->q;
  }
}
int Variable::parse_args(char *str, char **args)
{
  char *ptrnext;
  int narg = 0;
  char *ptr = str;
  while (ptr && narg < MAXFUNCARG) {
    ptrnext = find_next_comma(ptr);
    if (ptrnext) *ptrnext = '\0';
    args[narg] = utils::strdup(utils::trim(ptr));
    narg++;
    ptr = ptrnext;
    if (ptr) ptr++;
  }
  if (ptr) error->all(FLERR,"Too many args in variable function");
  return narg;
}
char *Variable::find_next_comma(char *str)
{
  int level = 0;
  for (char *p = str; *p; ++p) {
    if ('(' == *p) level++;
    else if (')' == *p) level--;
    else if (',' == *p && !level) return p;
  }
  return nullptr;
}
void Variable::print_var_error(const std::string &srcfile, const int lineno,
                               const std::string &errmsg, int ivar, int global)
{
  if ((ivar >= 0) && (ivar < nvar)) {
    std::string msg = fmt::format("Variable {}: ",names[ivar]) + errmsg;
    if (global)
      error->all(srcfile,lineno,msg);
    else
      error->one(srcfile,lineno,msg);
  } else {
    if (global)
      error->all(srcfile,lineno,errmsg);
    else
      error->one(srcfile,lineno,errmsg);
  }
}
