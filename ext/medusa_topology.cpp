#include "medusa_topology.h"
#include "jnc.h"

#include <deque>
#include <regex>

namespace medusa {

#define NEXT_WORD(msg)                                                                                                 \
  if (!next_word())                                                                                                    \
    throw std::runtime_error(error(msg));

#define NEXT_WORD_ASSERT(w, msg)                                                                                       \
  if (!next_word()) {                                                                                                  \
    throw std::runtime_error(error(msg));                                                                              \
  }                                                                                                                    \
  if (word != w) {                                                                                                     \
    throw std::runtime_error(error(msg));                                                                              \
  }

class TopologyParser {
public:
  enum { COMMENT_BEGIN, COMMENT_END, COMMENT_LINE, SPACE, RETURN, WORD, SYMBOL, FILE_END }; // char type
  enum { NORMAL, ERROR, END };

  Topology *top;
  TopResidue *pres = NULL;
  std::map<std::string, int> amap; // map atom type name to id

  int iline = 1;
  std::string filename;
  std::ifstream ifile;

  char c; // char
  int ct; // char type
  bool has_char; // if the char is in the queue

  std::string word;
  std::list<std::string> queue; // word queue

  std::string error_message;

  TopologyParser(Topology *top_) : top(top_) {}

  void clear() {
    if (pres != NULL) {
      delete pres;
      pres = NULL;
    }
    if (ifile.is_open()) {
      ifile.close();
    }
  }

  ~TopologyParser() { clear(); }

  std::string error(const std::string &msg) {
    throw std::runtime_error(jnc::string_format("File %s Line %d: %s.", filename, iline, msg));
  }

  void print_error_message() { std::cerr << error_message << std::endl; }

  int char_type(char c) {
    switch (c) {
    case '{':
      return COMMENT_BEGIN;
    case '}':
      return COMMENT_END;
    case '!':
      return COMMENT_LINE;
    case ' ':
      return SPACE;
    case '=':
      return SYMBOL;
    case '\n':
    case '\r':
      return RETURN;
    default:
      return WORD;
    }
  }

  void preview_char() {
    if (!has_char) {
      if (ifile) {
        ifile.get(c);
        ct = char_type(c);
        if (c == '\n')
          iline++;
      } else {
        c = '\0';
        ct = FILE_END;
      }
      has_char = true;
    }
  }

  void parse_spaces() {
    preview_char();
    if (ct == SPACE || ct == RETURN) {
      fetch_char();
      while (true) {
        preview_char();
        if (ct != SPACE && ct != RETURN) {
          return;
        } else {
          fetch_char();
        }
      }
    }
  }

  void parse_comment_m() {
    preview_char();
    if (ct == COMMENT_BEGIN) {
      fetch_char();
      int comment_level = 1;
      while (true) {
        fetch_char();
        if (ct == COMMENT_BEGIN) {
          comment_level++;
        } else if (ct == COMMENT_END) {
          comment_level--;
        } else if (ct == FILE_END) {
          error("The comment is not integral!");
          return;
        }

        if (comment_level == 0) {
          return;
        } else if (comment_level < 0) {
          error("The format of the comment is wrong");
          return;
        }
      }
    }
  }

  void parse_comment_s() {
    preview_char();
    if (ct == COMMENT_LINE) {
      fetch_char();
      while (true) {
        fetch_char();
        if (ct == RETURN || ct == FILE_END) {
          return;
        }
      }
    }
  }

  std::string parse_symbol() {
    preview_char();
    if (ct == SYMBOL) {
      std::stringstream stream;
      do {
        fetch_char();
        stream << c;
        preview_char();
      } while (ct == SYMBOL);
      return stream.str();
    }
  }

  std::string parse_word() {
    preview_char();
    if (ct == WORD) {
      std::stringstream stream;
      do {
        fetch_char();
        stream << c;
        preview_char();
      } while (ct == WORD);
      return stream.str();
    }
  }

  void fetch_char() {
    preview_char();
    has_char = false;
  }

  std::string next_word() {
    while (true) {
      preview_char();
      if (ct == SPACE || ct == RETURN) {
        parse_spaces();
      } else if (ct == COMMENT_BEGIN) {
        parse_comment_m();
      } else if (ct == COMMENT_END) {
        error("Wrong '}'");
      } else if (ct == COMMENT_LINE) {
        parse_comment_s();
      } else if (ct == SYMBOL) {
        return parse_symbol();
      } else if (ct == WORD) {
        return parse_word();
      } else if (ct == FILE_END) {
        return "";
      }
    }
  }

  std::string &preview() {
    if (queue.empty()) {
      queue.push_back(next_word());
    }
    return queue.front();
  }

  void fetch(const std::string &msg = "error") {
    preview();
    word = std::move(queue.front());
    if (word == "")
      error(msg);
    queue.pop_front();
  }

  void fetch_required(const std::string &required, const std::string &msg = "error") {
    fetch(msg);
    if (word != required)
      error(msg);
  }

  void parse_atom_type() {
    if (preview() == "MASS") {
      auto &atable = top->atom_table;
      fetch_required("MASS");
      fetch("Atom Type needs a name");
      auto &&atom_name = std::move(word);
      fetch("Atom Type needs a mass");
      double atom_mass = JN_DBL(word);
      if (std::none_of(atable.begin(), atable.end(),
                       [&atom_name](const TopAtomType &a) { return a.name == atom_name; })) {
        int id = atable.size();
        amap[atom_name] = id;
        atable.push_back(TopAtomType{.id = id, .mass = atom_mass, .name = atom_name});
      }
    }
  }

  void parse_atom_types() {
    while (true) {
      auto &n = preview();
      if (n == "MASS") {
        parse_atom_type();
      } else if (n == "END") {
        fetch();
        return;
      } else if (n == "") {
        return;
      }
    }
  }

  void parse_atom() {
    if (preview() == "ATOM") {
      fetch();

      fetch("Atom needs a name!");
      TopAtom atom;
      atom.name = std::move(word);

      fetch_required("TYPE", "Keyword 'TYPE' is needed");
      fetch_required("=", "'=' is needed");
      fetch("atom type is needed");
      auto &&atom_type = std::move(word);
      if (amap.find(atom_type) == amap.end()) {
        error(jnc::string_format("Atom type %s not defined", atom_type));
      } else {
        atom.type = amap[atom_type];
      }

      fetch_required("CHARge", "'CHARge' is needed");
      fetch_required("=", "'=' is needed");
      fetch("charge value is needed");
      atom.charge = JN_DBL(word);
      pres->atoms.push_back(std::move(atom));

      fetch_required("END", "'END' is needed");
    }
  }

  void parse_bond() {
    if (preview() == "BOND") {
      fetch();
      auto &v = pres->atoms;
      std::array<int, 2> bond;
      for (int i = 0; i < 2; i++) {
        fetch("bond needs two atoms");
        bond[i] = std::distance(v.begin(), std::find(v.begin(), v.end(), word));
      }
      pres->bonds.push_back(std::move(bond));
    }
  }

  void parse_dihedral() {
    if (preview() == "DIHE") {
      fetch();
      auto &v = pres->atoms;
      std::array<int, 4> dihedral;
      for (int i = 0; i < 4; i++) {
        fetch("dihedral needs four atoms");
        dihedral[i] = std::distance(v.begin(), std::find(v.begin(), v.end(), word));
      }
      pres->dihedrals.push_back(std::move(dihedral));
    }
  }

  void parse_improper() {
    if (preview() == "IMPROPER") {
      fetch();
      auto &v = pres->atoms;
      std::array<int, 4> improper;
      for (int i = 0; i < 4; i++) {
        fetch("improper needs four atoms");
        improper[i] = std::distance(v.begin(), std::find(v.begin(), v.end(), word));
      }
      pres->impropers.push_back(std::move(improper));
    }
  }

  void parse_residue() {
    if (preview() == "RESIdue") {
      fetch();
      fetch("Residue needs a name");
      pres = new TopResidue;
      pres->name = word;
      while (true) {
        auto &n = preview();
        if (n == "GROUP") {
          fetch();
        } else if (n == "ATOM") {
          parse_atom();
        } else if (n == "BOND") {
          parse_bond();
        } else if (n == "DIHE") {
          parse_dihedral();
        } else if (n == "END") {
          fetch();
          return;
        } else {
          error(jnc::string_format("'%s' may be wrong!", n));
        }
      }
    }
  }

  void parse_residues() {
    while (true) {
      auto &n = preview();
      if (n == "RESIdue") {
        parse_residue();
      } else if (n == "") {
        return;
      }
    }
  }

  int parse(const std::string &filename) {
    int r;
    try {
      ifile.open(filename.c_str());
      parse_atom_types();
      parse_residues();
      r = 1;
    } catch (const std::runtime_error &e) {
      std::cerr << e.what() << std::endl;
      r = 0;
    }
    clear();
    return r;
  }
};

void Topology::init(const std::string &name) {
  TopologyParser parser(this);
  parser.parse(name);
}

} // namespace medusa
