#ifndef LMP_TOKENIZER_H
#define LMP_TOKENIZER_H
#include "lmptype.h"
#include <exception>
#include <string>
#include <vector>
namespace LAMMPS_NS {
#define TOKENIZER_DEFAULT_SEPARATORS " \t\r\n\f"
class Tokenizer {
  std::string text;
  std::string separators;
  size_t start;
  size_t ntokens;

public:
  Tokenizer(std::string str,
            std::string separators = TOKENIZER_DEFAULT_SEPARATORS);
  Tokenizer(Tokenizer &&);
  Tokenizer(const Tokenizer &);
  Tokenizer &operator=(const Tokenizer &);
  Tokenizer &operator=(Tokenizer &&);
  void swap(Tokenizer &);
  void reset();
  void skip(int n = 1);
  bool has_next() const;
  bool contains(const std::string &str) const;
  std::string next();
  size_t count();
  std::vector<std::string> as_vector();
};
class TokenizerException : public std::exception {
  std::string message;

public:
  TokenizerException() = delete;
  explicit TokenizerException(const std::string &msg, const std::string &token);
  const char *what() const noexcept override { return message.c_str(); }
};
class InvalidIntegerException : public TokenizerException {
public:
  explicit InvalidIntegerException(const std::string &token)
      : TokenizerException("Not a valid integer number", token) {}
};
class InvalidFloatException : public TokenizerException {
public:
  explicit InvalidFloatException(const std::string &token)
      : TokenizerException("Not a valid floating-point number", token) {}
};
class ValueTokenizer {
  Tokenizer tokens;

public:
  ValueTokenizer(const std::string &str,
                 const std::string &separators = TOKENIZER_DEFAULT_SEPARATORS);
  ValueTokenizer(const ValueTokenizer &) = default;
  ValueTokenizer(ValueTokenizer &&);
  ValueTokenizer &operator=(const ValueTokenizer &);
  ValueTokenizer &operator=(ValueTokenizer &&);
  void swap(ValueTokenizer &);
  std::string next_string();
  tagint next_tagint();
  bigint next_bigint();
  int next_int();
  double next_double();
  bool has_next() const;
  bool contains(const std::string &value) const;
  void skip(int ntokens = 1);
  size_t count();
};
} // namespace LAMMPS_NS
#endif
