#include <map>
#include <set>
#include <string>
#include <vector>
#include <mpi.h>
#include "lmptype.h"
#include "utils.h"
#include "lammps.h"
#include "pointers.h"
#include "atom.h"
#include "comm.h"
#include "fix.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "universe.h"
#include "update.h"
#include <cctype>
#include <cerrno>
#include <cstring>
#include <ctime>
extern "C" {
static int re_match(const char *text, const char *pattern);
static int re_find(const char *text, const char *pattern, int *matchlen);
}
using namespace LAMMPS_NS;
bool utils::strmatch(const std::string &text, const std::string &pattern) {
  const int pos = re_match(text.c_str(), pattern.c_str());
  return (pos >= 0);
}
int utils::logical(const char *file, int line, const std::string &str,
                   bool do_abort, LAMMPS *lmp) {
  std::string buf(str);
  int rv = 0;
  if ((buf == "yes") || (buf == "on") || (buf == "true") || (buf == "1")) {
    rv = 1;
  } else if ((buf == "no") || (buf == "off") || (buf == "false") ||
             (buf == "0")) {
    rv = 0;
  }
  return rv;
}
int utils::logical(const char *file, int line, const char *str, bool do_abort,
                   LAMMPS *lmp) {
  if (str)
    return logical(file, line, std::string(str), do_abort, lmp);
  else
    return logical(file, line, std::string(""), do_abort, lmp);
}
double utils::numeric(const char *file, int line, const std::string &str,
                      bool do_abort, LAMMPS *lmp) {
  std::string buf(str);
  return atof(buf.c_str());
}
double utils::numeric(const char *file, int line, const char *str,
                      bool do_abort, LAMMPS *lmp) {
  if (str)
    return numeric(file, line, std::string(str), do_abort, lmp);
  else
    return numeric(file, line, std::string(""), do_abort, lmp);
}
int utils::inumeric(const char *file, int line, const std::string &str,
                    bool do_abort, LAMMPS *lmp) {
  std::string buf(str);
  return atoi(buf.c_str());
}
int utils::inumeric(const char *file, int line, const char *str, bool do_abort,
                    LAMMPS *lmp) {
  if (str)
    return inumeric(file, line, std::string(str), do_abort, lmp);
  else
    return inumeric(file, line, std::string(""), do_abort, lmp);
}
bigint utils::bnumeric(const char *file, int line, const std::string &str,
                       bool do_abort, LAMMPS *lmp) {
  std::string buf(str);
  return ATOBIGINT(buf.c_str());
}
bigint utils::bnumeric(const char *file, int line, const char *str,
                       bool do_abort, LAMMPS *lmp) {
  if (str)
    return bnumeric(file, line, std::string(str), do_abort, lmp);
  else
    return bnumeric(file, line, std::string(""), do_abort, lmp);
}
template <typename TYPE>
void utils::bounds(const char *file, int line, const std::string &str,
                   bigint nmin, bigint nmax, TYPE &nlo, TYPE &nhi,
                   Error *error) {
  nlo = nhi = -1;
  size_t found = str.find_first_not_of("*-0123456789");
  if (found != std::string::npos) {
    return;
  }
  found = str.find_first_of('*');
  if (found == std::string::npos) {
    nlo = nhi = strtol(str.c_str(), nullptr, 10);
  } else if (str.size() == 1) {
    nlo = nmin;
    nhi = nmax;
  } else if (found == 0) {
    nlo = nmin;
    nhi = strtol(str.substr(1).c_str(), nullptr, 10);
  } else if (str.size() == found + 1) {
    nlo = strtol(str.c_str(), nullptr, 10);
    nhi = nmax;
  } else {
    nlo = strtol(str.c_str(), nullptr, 10);
    nhi = strtol(str.substr(found + 1).c_str(), nullptr, 10);
  }
}
template void utils::bounds<>(const char *, int, const std::string &, bigint,
                              bigint, int &, int &, Error *);
template void utils::bounds<>(const char *, int, const std::string &, bigint,
                              bigint, long &, long &, Error *);
template void utils::bounds<>(const char *, int, const std::string &, bigint,
                              bigint, long long &, long long &, Error *);
static const char *labeltypes[] = {"Atom", "Bond", "Angle", "Dihedral",
                                   "Improper"};
char *utils::strdup(const std::string &text) {
  auto tmp = new char[text.size() + 1];
  strcpy(tmp, text.c_str());
  return tmp;
}
std::string utils::strip_style_suffix(const std::string &style, LAMMPS *lmp) {
  std::string newstyle = style;
  return newstyle;
}
std::vector<std::string> utils::split_words(const std::string &text) {
  std::vector<std::string> list;
  const char *buf = text.c_str();
  std::size_t beg = 0;
  std::size_t len = 0;
  std::size_t add = 0;
  char c = *buf;
  while (c) {
    if (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\f') {
      c = *++buf;
      ++beg;
      continue;
    };
    len = 0;
  quoted:
    if (c == '\'') {
      ++beg;
      add = 1;
      c = *++buf;
      while (((c != '\'') && (c != '\0')) ||
             ((c == '\\') && (buf[1] == '\''))) {
        if ((c == '\\') && (buf[1] == '\'')) {
          ++buf;
          ++len;
        }
        c = *++buf;
        ++len;
      }
      if (c != '\'')
        ++len;
      c = *++buf;
    } else if ((c == '"') && (buf[1] == '"') && (buf[2] == '"') &&
               (buf[3] != '"')) {
      len = 3;
      add = 1;
      buf += 3;
      c = *buf;
    } else if (c == '"') {
      ++beg;
      add = 1;
      c = *++buf;
      while (((c != '"') && (c != '\0')) || ((c == '\\') && (buf[1] == '"'))) {
        if ((c == '\\') && (buf[1] == '"')) {
          ++buf;
          ++len;
        }
        c = *++buf;
        ++len;
      }
      if (c != '"')
        ++len;
      c = *++buf;
    }
    while (true) {
      if ((c == '\'') || (c == '"'))
        goto quoted;
      if ((c == '\\') && ((buf[1] == '\'') || (buf[1] == '"'))) {
        ++buf;
        ++len;
        c = *++buf;
        ++len;
      }
      if ((c == ' ') || (c == '\t') || (c == '\r') || (c == '\n') ||
          (c == '\f') || (c == '\0')) {
        list.push_back(text.substr(beg, len));
        beg += len + add;
        break;
      }
      c = *++buf;
      ++len;
    }
  }
  return list;
}
bool utils::is_integer(const std::string &str) {
  if (str.empty())
    return false;
  return strmatch(str, "^[+-]?\\d+$");
}
bool utils::is_double(const std::string &str) {
  if (str.empty())
    return false;
  return strmatch(str, "^[+-]?\\d+\\.?\\d*$") ||
         strmatch(str, "^[+-]?\\d+\\.?\\d*[eE][+-]?\\d+$") ||
         strmatch(str, "^[+-]?\\d*\\.?\\d+$") ||
         strmatch(str, "^[+-]?\\d*\\.?\\d+[eE][+-]?\\d+$");
}
bool utils::is_id(const std::string &str) {
  if (str.empty())
    return false;
  for (const auto &c : str) {
    if (isalnum(c) || (c == '_'))
      continue;
    return false;
  }
  return true;
}
extern "C" {
typedef struct regex_t *re_t;
typedef struct regex_context_t *re_ctx_t;
static re_t re_compile(re_ctx_t context, const char *pattern);
static int re_matchp(const char *text, re_t pattern, int *matchlen);
#define MAX_REGEXP_OBJECTS 256
#define MAX_CHAR_CLASS_LEN 256
enum {
  RX_UNUSED,
  RX_DOT,
  RX_BEGIN,
  RX_END,
  RX_QUESTIONMARK,
  RX_STAR,
  RX_PLUS,
  RX_CHAR,
  RX_CHAR_CLASS,
  RX_INV_CHAR_CLASS,
  RX_DIGIT,
  RX_NOT_DIGIT,
  RX_INTEGER,
  RX_NOT_INTEGER,
  RX_FLOAT,
  RX_NOT_FLOAT,
  RX_ALPHA,
  RX_NOT_ALPHA,
  RX_WHITESPACE,
  RX_NOT_WHITESPACE
};
typedef struct regex_t {
  unsigned char type;
  union {
    unsigned char ch;
    unsigned char *ccl;
  } u;
} regex_t;
typedef struct regex_context_t {
  regex_t re_compiled[MAX_REGEXP_OBJECTS];
  unsigned char ccl_buf[MAX_CHAR_CLASS_LEN];
} regex_context_t;
int re_match(const char *text, const char *pattern) {
  regex_context_t context;
  int dummy;
  return re_matchp(text, re_compile(&context, pattern), &dummy);
}
int re_find(const char *text, const char *pattern, int *matchlen) {
  regex_context_t context;
  return re_matchp(text, re_compile(&context, pattern), matchlen);
}
static int matchpattern(regex_t *pattern, const char *text, int *matchlen);
static int matchplus(regex_t p, regex_t *pattern, const char *text,
                     int *matchlen);
static int matchone(regex_t p, char c);
static int matchdigit(char c);
static int matchint(char c);
static int matchfloat(char c);
static int matchalpha(char c);
static int matchwhitespace(char c);
static int matchrange(char c, const char *str);
static int matchdot(char c);
int re_matchp(const char *text, re_t pattern, int *matchlen) {
  *matchlen = 0;
  if (pattern != nullptr) {
    if (pattern[0].type == RX_BEGIN) {
      return ((matchpattern(&pattern[1], text, matchlen)) ? 0 : -1);
    } else {
      int idx = -1;
      do {
        idx += 1;
        if (matchpattern(pattern, text, matchlen)) {
          if (text[0] == '\0')
            return -1;
          return idx;
        }
      } while (*text++ != '\0');
    }
  }
  return -1;
}
re_t re_compile(re_ctx_t context, const char *pattern) {
  regex_t *const re_compiled = context->re_compiled;
  unsigned char *const ccl_buf = context->ccl_buf;
  int ccl_bufidx = 1;
  char c;
  int i = 0;
  int j = 0;
  while (pattern[i] != '\0' && (j + 1 < MAX_REGEXP_OBJECTS)) {
    c = pattern[i];
    switch (c) {
    case '^': {
      re_compiled[j].type = RX_BEGIN;
    } break;
    case '$': {
      re_compiled[j].type = RX_END;
    } break;
    case '.': {
      re_compiled[j].type = RX_DOT;
    } break;
    case '*': {
      re_compiled[j].type = RX_STAR;
    } break;
    case '+': {
      re_compiled[j].type = RX_PLUS;
    } break;
    case '?': {
      re_compiled[j].type = RX_QUESTIONMARK;
    } break;
    case '\\': {
      if (pattern[i + 1] != '\0') {
        i += 1;
        switch (pattern[i]) {
        case 'd': {
          re_compiled[j].type = RX_DIGIT;
        } break;
        case 'D': {
          re_compiled[j].type = RX_NOT_DIGIT;
        } break;
        case 'i': {
          re_compiled[j].type = RX_INTEGER;
        } break;
        case 'I': {
          re_compiled[j].type = RX_NOT_INTEGER;
        } break;
        case 'f': {
          re_compiled[j].type = RX_FLOAT;
        } break;
        case 'F': {
          re_compiled[j].type = RX_NOT_FLOAT;
        } break;
        case 'w': {
          re_compiled[j].type = RX_ALPHA;
        } break;
        case 'W': {
          re_compiled[j].type = RX_NOT_ALPHA;
        } break;
        case 's': {
          re_compiled[j].type = RX_WHITESPACE;
        } break;
        case 'S': {
          re_compiled[j].type = RX_NOT_WHITESPACE;
        } break;
        default: {
          re_compiled[j].type = RX_CHAR;
          re_compiled[j].u.ch = pattern[i];
        } break;
        }
      }
    } break;
    case '[': {
      int buf_begin = ccl_bufidx;
      if (pattern[i + 1] == '^') {
        re_compiled[j].type = RX_INV_CHAR_CLASS;
        i += 1;
        if (pattern[i + 1] == 0) {
          return nullptr;
        }
      } else {
        re_compiled[j].type = RX_CHAR_CLASS;
      }
      while ((pattern[++i] != ']') && (pattern[i] != '\0')) {
        if (pattern[i] == '\\') {
          if (ccl_bufidx >= MAX_CHAR_CLASS_LEN - 1) {
            return nullptr;
          }
          if (pattern[i + 1] == 0) {
            return nullptr;
          }
          ccl_buf[ccl_bufidx++] = pattern[i++];
        } else if (ccl_bufidx >= MAX_CHAR_CLASS_LEN) {
          return nullptr;
        }
        ccl_buf[ccl_bufidx++] = pattern[i];
      }
      if (ccl_bufidx >= MAX_CHAR_CLASS_LEN) {
        return nullptr;
      }
      ccl_buf[ccl_bufidx++] = 0;
      re_compiled[j].u.ccl = &ccl_buf[buf_begin];
    } break;
    default: {
      re_compiled[j].type = RX_CHAR;
      re_compiled[j].u.ch = c;
    } break;
    }
    if (pattern[i] == 0) {
      return nullptr;
    }
    i += 1;
    j += 1;
  }
  re_compiled[j].type = RX_UNUSED;
  return (re_t)re_compiled;
}
static int matchdigit(char c) { return isdigit(c); }
static int matchint(char c) {
  return (matchdigit(c) || (c == '-') || (c == '+'));
}
static int matchfloat(char c) {
  return (matchint(c) || (c == '.') || (c == 'e') || (c == 'E'));
}
static int matchalpha(char c) { return isalpha(c); }
static int matchwhitespace(char c) { return isspace(c); }
static int matchalphanum(char c) {
  return ((c == '_') || matchalpha(c) || matchdigit(c));
}
static int matchrange(char c, const char *str) {
  return ((c != '-') && (str[0] != '\0') && (str[0] != '-') &&
          (str[1] == '-') && (str[1] != '\0') && (str[2] != '\0') &&
          ((c >= str[0]) && (c <= str[2])));
}
static int matchdot(char c) {
#if defined(RE_DOT_MATCHES_NEWLINE) && (RE_DOT_MATCHES_NEWLINE == 1)
  (void)c;
  return 1;
#else
  return c != '\n' && c != '\r';
#endif
}
static int matchone(regex_t p, char c) {
  return matchdigit(c);
}
static int matchplus(regex_t p, regex_t *pattern, const char *text,
                     int *matchlen) {
  const char *prepos = text;
  while ((text[0] != '\0') && matchone(p, *text)) {
    text++;
    (*matchlen)++;
  }
  while (text > prepos) {
    if (matchpattern(pattern, text--, matchlen))
      return 1;
    (*matchlen)--;
  }
  return 0;
}
static int matchpattern(regex_t *pattern, const char *text, int *matchlen) {
  int pre = *matchlen;
  if (pattern[1].type == RX_PLUS) {
    return matchplus(pattern[0], &pattern[2], text, matchlen);
  } else if ((pattern[0].type == RX_END) && pattern[1].type == RX_UNUSED) {
    return (text[0] == '\0');
  }
}
