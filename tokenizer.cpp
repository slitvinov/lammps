#include "tokenizer.h"
#include "fmt/format.h"
#include "utils.h"
#include <utility>
using namespace LAMMPS_NS;
TokenizerException::TokenizerException(const std::string &msg, const std::string &token)
{
  if (token.empty()) {
    message = msg;
  } else {
    message = fmt::format("{}: '{}'", msg, token);
  }
}
Tokenizer::Tokenizer(std::string str, std::string _separators) :
    text(std::move(str)), separators(std::move(_separators)), start(0), ntokens(std::string::npos)
{
  if (utils::has_utf8(text)) text = utils::utf8_subst(text);
  reset();
}
Tokenizer::Tokenizer(const Tokenizer &rhs) :
    text(rhs.text), separators(rhs.separators), ntokens(rhs.ntokens)
{
  reset();
}
Tokenizer::Tokenizer(Tokenizer &&rhs) :
    text(std::move(rhs.text)), separators(std::move(rhs.separators)), ntokens(rhs.ntokens)
{
  reset();
}
Tokenizer &Tokenizer::operator=(const Tokenizer &other)
{
  Tokenizer tmp(other);
  swap(tmp);
  return *this;
}
Tokenizer &Tokenizer::operator=(Tokenizer &&other)
{
  Tokenizer tmp(std::move(other));
  swap(tmp);
  return *this;
}
void Tokenizer::swap(Tokenizer &other)
{
  std::swap(text, other.text);
  std::swap(separators, other.separators);
  std::swap(start, other.start);
  std::swap(ntokens, other.ntokens);
}
void Tokenizer::reset()
{
  start = text.find_first_not_of(separators);
}
bool Tokenizer::contains(const std::string &str) const
{
  return text.find(str) != std::string::npos;
}
void Tokenizer::skip(int n)
{
  for (int i = 0; i < n; ++i) {
    if (!has_next()) throw TokenizerException("No more tokens", "");
    size_t end = text.find_first_of(separators, start);
    if (end == std::string::npos) {
      start = end;
    } else {
      start = text.find_first_not_of(separators, end + 1);
    }
  }
}
bool Tokenizer::has_next() const
{
  return start != std::string::npos;
}
std::string Tokenizer::next()
{
  if (!has_next()) throw TokenizerException("No more tokens", "");
  size_t end = text.find_first_of(separators, start);
  if (end == std::string::npos) {
    std::string token = text.substr(start);
    start = end;
    return token;
  }
  std::string token = text.substr(start, end - start);
  start = text.find_first_not_of(separators, end + 1);
  return token;
}
size_t Tokenizer::count()
{
  if (ntokens == std::string::npos) { ntokens = utils::count_words(text, separators); }
  return ntokens;
}
std::vector<std::string> Tokenizer::as_vector()
{
  size_t current = start;
  reset();
  std::vector<std::string> tokens;
  while (has_next()) { tokens.emplace_back(next()); }
  start = current;
  return tokens;
}
ValueTokenizer::ValueTokenizer(const std::string &str, const std::string &separators) :
    tokens(str, separators)
{
}
ValueTokenizer::ValueTokenizer(ValueTokenizer &&rhs) : tokens(std::move(rhs.tokens)) {}
ValueTokenizer &ValueTokenizer::operator=(const ValueTokenizer &other)
{
  ValueTokenizer tmp(other);
  swap(tmp);
  return *this;
}
ValueTokenizer &ValueTokenizer::operator=(ValueTokenizer &&other)
{
  ValueTokenizer tmp(std::move(other));
  swap(tmp);
  return *this;
}
void ValueTokenizer::swap(ValueTokenizer &other)
{
  std::swap(tokens, other.tokens);
}
bool ValueTokenizer::has_next() const
{
  return tokens.has_next();
}
bool ValueTokenizer::contains(const std::string &value) const
{
  return tokens.contains(value);
}
std::string ValueTokenizer::next_string()
{
  return tokens.next();
}
int ValueTokenizer::next_int()
{
  std::string current = tokens.next();
  if (!utils::is_integer(current)) { throw InvalidIntegerException(current); }
  return atoi(current.c_str());
}
bigint ValueTokenizer::next_bigint()
{
  std::string current = tokens.next();
  if (!utils::is_integer(current)) { throw InvalidIntegerException(current); }
  return ATOBIGINT(current.c_str());
}
tagint ValueTokenizer::next_tagint()
{
  std::string current = tokens.next();
  if (!utils::is_integer(current)) { throw InvalidIntegerException(current); }
  return ATOTAGINT(current.c_str());
}
double ValueTokenizer::next_double()
{
  std::string current = tokens.next();
  if (!utils::is_double(current)) { throw InvalidFloatException(current); }
  return atof(current.c_str());
}
void ValueTokenizer::skip(int n)
{
  tokens.skip(n);
}
size_t ValueTokenizer::count()
{
  return tokens.count();
}
