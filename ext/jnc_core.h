#pragma once

#include <cfloat>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <algorithm>
#include <functional>
#include <memory>
#include <type_traits>
#include <utility>

#include <array>
#include <deque>
#include <list>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "./Eigen/Core"
#include "./Eigen/Dense"
#include "./Eigen/Geometry"
#include "./Eigen/SVD"

namespace jnc {

constexpr double pi = 3.14159265358979323846;

template <typename T> using Vec = std::vector<T>;

enum { JN_CONTINUE, JN_BREAK, JN_GO };

template <typename Func_> void file_each_line(std::string fn, Func_ &&f) {
  std::ifstream ifile(fn.c_str());
  std::string line;
  while (ifile) {
    std::getline(ifile, line);
    auto r = f(line);
    if (r == JN_CONTINUE)
      continue;
    else if (r == JN_BREAK)
      break;
  }
  ifile.close();
}

#define JN_FILE_BEGIN(filename, line) \
{ \
  std::ifstream ifile((filename)); \
  std::string line; \
  while (ifile) { \
    std::getline(ifile, line); \

#define JN_FILE_END \
  } \
  ifile.close(); \
}

template <typename T, typename U> T lexical_cast(U &&u) {
  std::stringstream stream;
  T t;

  stream << u;
  stream >> t;
  return t;
}

#define JN_INT(a) jnc::lexical_cast<int>(a)
#define JN_DBL(a) jnc::lexical_cast<double>(a)
#define JN_FLT(a) jnc::lexical_cast<float>(a)
#define JN_STR(a) jnc::lexical_cast<std::string>(a)

template <typename T> void hash_combine_(std::size_t &seed, const T &value) {
  seed ^= (std::size_t)value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

template <class T> inline void hash_combine(std::size_t &seed, const T &v) {
  std::hash<T> hasher;
  hash_combine_(seed, hasher(v));
}

template <class It> std::size_t hash_range(It first, It last) {
  std::size_t seed = 0;
  for (; first != last; ++first) {
    hash_combine_(seed, *first);
  }
  return seed;
}

template <typename Stream_> void stream_push(Stream_ &&stream) {}

template <typename Stream_, typename Value_, typename... Values_>
void stream_push(Stream_ &&stream, Value_ &&value, Values_ &&...values) {
  stream << value;
  stream_push(stream, values...);
}

inline Vec<std::string> string_tokenize(const std::string &str, const std::string &delimiters) {
  Vec<std::string> tokens;
  auto lastPos = str.find_first_not_of(delimiters, 0);
  auto pos = str.find_first_of(delimiters, lastPos);
  while (std::string::npos != pos || std::string::npos != lastPos) {
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    lastPos = str.find_first_not_of(delimiters, pos);
    pos = str.find_first_of(delimiters, lastPos);
  }
  return tokens;
}

inline Vec<std::string> string_tokenize(const std::string &str, const std::string &delimiters,
                                const std::string &temp) {
  Vec<std::string> tokens;
  using pair_t = std::pair<std::string::size_type, std::string::size_type>;
  Vec<pair_t> vec;
  std::string::size_type first_i, first_j, second_i, second_j;
  int expected = 0;
  for (std::string::size_type i = 0; i < str.size(); i++) {
    int flag = 0;
    std::string::size_type j;
    for (j = 0; j < temp.size(); j++) {
      if (str[i] == temp[j]) {
        if (j % 2 == 0 && expected == 0) {
          flag = 1;
          break;
        } else if (j % 2 == 1 && expected == 1) {
          flag = 2;
          break;
        }
      }
    }
    if (flag == 1) {
      first_i = i;
      first_j = j;
      expected = 1;
    } else if (flag == 2 && j - first_j == 1) {
      second_i = i;
      second_j = j;
      expected = 0;
      vec.push_back(std::make_pair(first_i, second_i));
    }
  }
  auto lastPos = str.find_first_not_of(delimiters, 0);
  auto pos = str.find_first_of(delimiters, lastPos);
  while (std::any_of(vec.begin(), vec.end(), [&pos](const pair_t &p) {
    return pos != std::string::npos && p.first < pos && pos < p.second;
  })) {
    pos = str.find_first_of(delimiters, pos + 1);
  }
  while (std::string::npos != pos || std::string::npos != lastPos) {
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    lastPos = str.find_first_not_of(delimiters, pos);
    pos = str.find_first_of(delimiters, lastPos);
    while (std::any_of(vec.begin(), vec.end(), [&pos](const pair_t &p) {
      return pos != std::string::npos && p.first < pos && pos < p.second;
    })) {
      pos = str.find_first_of(delimiters, pos + 1);
    }
  }
  return tokens;
}

inline void string_upper(std::string &s) {
  std::transform(s.begin(), s.end(), s.begin(), ::toupper);
}

inline std::string string_upper_c(const std::string &str) {
  std::string s = str;
  std::transform(s.begin(), s.end(), s.begin(), ::toupper);
  return s;
}

inline void string_lower(std::string &s) {
  std::transform(s.begin(), s.end(), s.begin(), ::tolower);
}

inline std::string string_lower_c(const std::string &str) {
  std::string s = str;
  std::transform(s.begin(), s.end(), s.begin(), ::tolower);
  return s;
}

inline bool string_includes(const std::string &s1, const std::string &s2) {
  return s1.find(s2) != std::string::npos;
}

inline void string_trim(std::string &s) {
  if (!(s.empty())) {
    s.erase(0, s.find_first_not_of(" "));
    s.erase(s.find_last_not_of(" ") + 1);
  }
}

inline std::string string_trim_c(const std::string &s) {
  if (s.empty()) {
    return s;
  }
  std::string::size_type beg = s.find_first_not_of(" \t\n\r");
  std::string::size_type end = s.find_last_not_of(" \t\n\r");
  if (beg == std::string::npos || end == std::string::npos) {
    return "";
  } else {
    return std::string(s.begin() + beg, s.begin() + end + 1);
  }
}

inline bool string_starts_with(const std::string &str1, const std::string &str2) {
  int sz = str2.size();
  for (int i = 0; i < sz; i++) {
    if (str1[i] != str2[i])
      return false;
  }
  return true;
}

inline bool string_istarts_with(const std::string &str1, const std::string &str2) {
  int sz = str2.size();
  for (int i = 0; i < sz; i++) {
    if (tolower(str1[i]) != tolower(str2[i]))
      return false;
  }
  return true;
}

inline bool string_ends_with(const std::string &str1, const std::string &str2) {
  int sz = str2.size();
  auto it1 = str1.rbegin();
  auto it2 = str2.rbegin();
  for (int i = 0; i < sz; i++) {
    if (*it1 != *it2)
      return false;
    it1++;
    it2++;
  }
  return true;
}

inline bool string_iends_with(const std::string &str1, const std::string &str2) {
  int sz = str2.size();
  auto it1 = str1.rbegin();
  auto it2 = str2.rbegin();
  for (int i = 0; i < sz; i++) {
    if (tolower(*it1) != tolower(*it2))
      return false;
    it1++;
    it2++;
  }
  return true;
}

template <typename... Strs_> std::string string_merge(Strs_ &&...strs) {
  std::stringstream stream;
  stream_push(stream, strs...);
  return stream.str();
}

///////////////// string_format //////////////////
template <typename T> T string_to_chars(T &&t) { return t; }

inline const char *string_to_chars(std::string &s) { return s.c_str(); }

inline const char *string_to_chars(const std::string &s) { return s.c_str(); }

inline const char *string_to_chars(std::string &&s) { return s.c_str(); }

template <typename... Pars_> std::string string_format(std::string fmt, Pars_ &&...pars) {
  int count = snprintf(NULL, 0, fmt.c_str(),
                       string_to_chars(std::forward<Pars_>(pars))...);
  std::string buf;
  buf.resize(count);
  sprintf(&(buf[0]), fmt.c_str(), string_to_chars(pars)...);
  return buf;
}
//////////////////////////////////////////////////

template <typename _Interval, typename _Ls>
static std::string string_join(_Interval &&interval, _Ls &&ls) {
  std::stringstream stream;
  int i = 0;
  for (const auto &s : ls) {
    if (i != 0)
      stream << interval;
    stream << s;
    i++;
  }
  return stream.str();
}

inline bool string_iequals(const std::string &a, const std::string &b) {
  unsigned int sz = a.size();
  if (b.size() != sz)
    return false;
  for (unsigned int i = 0; i < sz; ++i)
    if (tolower(a[i]) != tolower(b[i]))
      return false;
  return true;
}

#define JN_DIE(...)                                                            \
  do {                                                                         \
    std::cerr << ::jnc::string_merge(__FILE__, "(", __LINE__,                  \
                                     "): ", __VA_ARGS__)                       \
              << std::endl;                                                    \
    std::abort();                                                              \
  } while (0)

#define JN_ASSERT(condition, what)                                             \
  if (!condition) {                                                            \
    JN_DIE(what);                                                              \
  }

template <typename V, int N> class ArrayHash {
public:
  std::size_t operator()(const std::array<V, N> &arr) const {
    std::size_t seed = 0;
    for (auto &&i : arr) {
      hash_combine(seed, i);
    }
    return seed;
  }
};


} // namespace jnc
