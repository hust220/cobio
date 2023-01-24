#pragma once

#include "./jnc_core.h"

namespace jnc {
namespace bio {

class PdbAtom : public std::array<double, 3> {
public:
  std::string name;
  std::string type;
  std::string element;
  int num;
  double charge;
  double bfactor = 0;
  bool is_std = true;
};

const static std::map<std::string, char> g_aa3 = {
    {"ALA", 'A'}, {"ARG", 'R'}, {"ASN", 'N'}, {"ASP", 'D'}, {"CYS", 'C'},
    {"GLN", 'Q'}, {"GLU", 'E'}, {"GLY", 'G'}, {"HIS", 'H'}, {"ILE", 'I'},
    {"LEU", 'L'}, {"LYS", 'K'}, {"MET", 'M'}, {"PHE", 'F'}, {"PRO", 'P'},
    {"SER", 'S'}, {"THR", 'T'}, {"TRP", 'W'}, {"TYR", 'Y'}, {"VAL", 'V'},
    {"PYL", 'O'}, {"SEC", 'U'}, {"ASX", 'B'}, {"GLX", 'Z'}, {"XAA", 'X'},
    {"XLE", 'J'}, {"TERM", '*'}};

class PdbResidue : public Vec<PdbAtom> {
public:
  std::string name;
  int num = -1;
  bool is_std = true;

  PdbAtom &atom(std::string name) {
    for (auto &&atom : *this) {
      if (atom.name == name)
        return atom;
    }
    JN_DIE("Atom '" + name + "' not found");
  }

  const PdbAtom &atom(std::string name) const {
    for (auto &&atom : *this) {
      if (atom.name == name)
        return atom;
    }
    JN_DIE("Atom '" + name + "' not found");
  }

  bool has_atom(std::string name) const {
    for (auto &&atom : *this) {
      if (atom.name == name)
        return true;
    }
    return false;
  }
};

class PdbChain : public Vec<PdbResidue> {
public:
  std::string name;

  PdbChain() {}

  Vec<const PdbAtom *> patoms() const {
    Vec<const PdbAtom *> as;
    for (auto &&res : *this) {
      for (auto &&atom : res) {
        as.push_back(&atom);
      }
    }
    return as;
  }

  Vec<PdbAtom *> patoms() {
    Vec<PdbAtom *> as;
    for (auto &&res : *this) {
      for (auto &&atom : res) {
        as.push_back(&atom);
      }
    }
    return as;
  }
};

class PdbModel : public Vec<PdbChain> {
public:
  std::string name;
  int num;

  Vec<const PdbAtom *> patoms() const {
    Vec<const PdbAtom *> as;
    for (auto &&chain : *this) {
      for (auto &&res : chain) {
        for (auto &&atom : res) {
          as.push_back(&atom);
        }
      }
    }
    return as;
  }

  Vec<PdbAtom *> patoms() {
    Vec<PdbAtom *> as;
    for (auto &&chain : *this) {
      for (auto &&res : chain) {
        for (auto &&atom : res) {
          as.push_back(&atom);
        }
      }
    }
    return as;
  }

  Vec<const PdbResidue *> presidues() const {
    Vec<const PdbResidue *> rs;
    for (auto &&chain : *this) {
      for (auto &&res : chain) {
        rs.push_back(&res);
      }
    }
    return std::move(rs);
  }

  Vec<PdbResidue *> presidues() {
    Vec<PdbResidue *> rs;
    int ir = 0;
    for (auto &&chain : *this) {
      for (auto &&res : chain) {
        // if (ir == 0) {
        //   std::cout << "presidues" << std::endl;
        //   std::cout << &res << std::endl;
        // }
        rs.push_back(&res);
        ir++;
      }
    }
    return std::move(rs);
  }
};

class Pdb;
Pdb read_pdb(const std::string &filename, Pdb &);

class Pdb : public Vec<PdbModel> {
public:
  std::string name;

  Pdb() {}

  Pdb(const std::string &filename) { read_pdb(filename, *this); }

  void remove_hydrogens() {
    for (auto &&model : *this) {
      for (auto &&chain : model) {
        for (auto &&res : chain) {
          auto r = res;
          r.clear();
          for (auto &&atom : res) {
            //                    if (atom.name[0] != 'H' &&
            //                    (!std::isdigit(atom.name[0]) || atom.name[1]
            //                    != 'H'))
            if (atom.name[0] != 'H' && atom.name[1] != 'H') {
              r.push_back(atom);
            }
          }
          res = std::move(r);
        }
      }
    }
  }

  void sort() {
    for (auto &&model : *this) {
      for (auto &&chain : model) {
        for (auto &&res : chain) {
          std::sort(res.begin(), res.end(), [](const PdbAtom &a1, const PdbAtom &a2) {
            return a1.name < a2.name;
          });
        }
      }
    }
  }
};

/**
 * Write a line in a PDB file.
 */
inline void write_line(std::ostream &out, double x, double y, double z,
                       const std::string &atom_name, int atom_index,
                       const std::string &residue_name, int residue_index,
                       const std::string &chain_name) {
  out << string_format(
             "ATOM%7d  %-4.4s%3.3s%2.2s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%12.12s  ",
             atom_index + 1, atom_name, residue_name, chain_name,
             residue_index + 1, x, y, z, 1.00, 0.00, "" + std::string(1, atom_name[0]))
      << std::endl;
}

class PdbReader {
public:
  PdbResidue atoms;
  PdbChain residues;
  PdbModel chains;
  Pdb &models;

  PdbReader(Pdb &pdb) : models(pdb) {}

  int model_num = 0;

  bool continue_reading = true;

  struct ParsedLine {
    std::string atom_name, atom_type, atom_flag, atom_element, res_name, res_flag,
        chain_name;
    int atom_num, res_num;
    double x, y, z, atom_bfactor;
    bool is_std;
  };

  ParsedLine ol;

  void parse_line(std::string line, ParsedLine &pl) {
    pl.res_name = string_trim_c(line.substr(17, 3));
    pl.res_num = JN_INT(string_trim_c(line.substr(22, 4)));
    pl.res_flag = string_trim_c(line.substr(26, 1));
    pl.chain_name = string_trim_c(line.substr(20, 2));
    pl.is_std = (!line.compare(0, 4, "ATOM"));

    pl.atom_name = string_trim_c(line.substr(12, 4));
    pl.atom_num = JN_INT(string_trim_c(line.substr(6, 5)));
    pl.atom_flag = string_trim_c(line.substr(16, 1));

    if (line.size() >= 66) {
      pl.atom_bfactor = JN_DBL(line.substr(60, 6));
    } else {
      pl.atom_bfactor = 0.0;
    }

    if (line.size() >= 78) {
      pl.atom_element = line.substr(77, 1);
    } else {
      pl.atom_element = "X";
    }

    pl.x = JN_DBL(string_trim_c(line.substr(30, 8)));
    pl.y = JN_DBL(string_trim_c(line.substr(38, 8)));
    pl.z = JN_DBL(string_trim_c(line.substr(46, 8)));

    if (line.size() >= 78) {
      pl.atom_type = string_trim_c(line.substr(76, 2));
    } else {
      pl.atom_type = "X";
    }

    std::replace(pl.atom_name.begin(), pl.atom_name.end(), '*', '\'');
    if (pl.atom_name == "O1P")
      pl.atom_name = "OP1";
    if (pl.atom_name == "O2P")
      pl.atom_name = "OP2";
  };

  void add_residue() {
    if (!atoms.empty()) {
      PdbResidue residue = std::move(atoms);
      residue.is_std =
          std::all_of(atoms.begin(), atoms.end(),
                      [](const PdbAtom &atom) { return atom.is_std; });
      residue.name = ol.res_name;
      residue.num = ol.res_num;
      residues.push_back(std::move(residue));
    }
  };

  void add_chain() {
    add_residue();
    if (!residues.empty()) {
      PdbChain chain = std::move(residues);
      chain.name = ol.chain_name;
      chains.push_back(std::move(chain));
    }
  };

  void add_atom(std::string line) {
    ParsedLine pl;
    parse_line(line, pl);

    if (!atoms.empty()) {
      if (pl.res_num != ol.res_num || pl.res_name != ol.res_name ||
          pl.res_flag != ol.res_flag || pl.chain_name != ol.chain_name) {
        add_residue();
      }
      if (pl.chain_name != ol.chain_name) {
        add_chain();
      }
    }

    //        if (string_starts_with(pl.atom_name, "H")) return;
    if (pl.res_name == "HOH" || pl.res_name == "H2O")
      return;

    // Check whether the atom has appeared before.
    if (!pl.atom_flag.empty()) {
      if (std::any_of(atoms.begin(), atoms.end(), [&pl](const PdbAtom &atom) {
            return atom.name == pl.atom_name;
          }))
        return;
    }

    PdbAtom atom;
    atom[0] = pl.x;
    atom[1] = pl.y;
    atom[2] = pl.z;
    atom.name = pl.atom_name;
    atom.type = pl.atom_type;
    atom.num = pl.atom_num;
    atom.element = pl.atom_element;
    atom.is_std = pl.is_std;
    atom.bfactor = pl.atom_bfactor;

    atoms.push_back(std::move(atom));

    ol = pl;
  };

  void add_model() {
    add_chain();
    if (!chains.empty()) {
      PdbModel model = std::move(chains);
      model.num = model_num;
      models.push_back(std::move(model));
      model_num++;
    }
  };

  void read_line(const std::string &line) {
    if (continue_reading) {
      if (string_starts_with(line, "ATOM")) {
        add_atom(line);
      } else if (string_starts_with(line, "HETATM")) {
        add_atom(line);
      } else if (string_starts_with(line, "MODEL")) {
        add_model();

        auto &&v = string_tokenize(line, " ");
        if (v.size() >= 2) {
          model_num = JN_INT(v[1]) - 1;
        }
      } else if (string_starts_with(line, "ENDMDL")) {
        add_model();
      } else if (string_starts_with(line, "TER")) {
        add_chain();
      } else if (string_trim_c(line) == "END") {
        continue_reading = false;
      } else {
        //
      }
    }
  }

  void end_reading() {
    add_model();
    continue_reading = false;
  }

  void read(const std::string &fn) {
    std::ifstream ifile(fn.c_str());
    continue_reading = true;
    while (ifile) {
      std::string line;
      std::getline(ifile, line);
      read_line(line);
    }
    ifile.close();

    end_reading();
    models.name = fn;
  }
};

inline Pdb read_pdb(const std::string &filename, Pdb &pdb) {
  PdbReader pdb_reader(pdb);
  pdb_reader.read(filename);
  return pdb;
}

inline Pdb read_pdb(const std::string &filename) {
  Pdb pdb;
  return read_pdb(filename, pdb);
}

class PdbWriter {
public:
  const Vec<std::string> chain_names = {"A", "B", "C", "D", "E", "F", "G", "H", "I",
                                "J", "K", "L", "M", "N", "O", "P", "Q", "R",
                                "S", "T", "U", "V", "W", "X", "Y", "Z"};

  int atom_num, residue_num, model_num;
  std::string atom_name, residue_name, chain_name, atom_name_label;
  double x, y, z, a, b;

  std::ostream stream;

  PdbWriter() : stream(std::cout.rdbuf()) { init(); }

  PdbWriter(std::ostream &out) : stream(out.rdbuf()) { init(); }

  void init() {
    atom_num = 1;
    residue_num = 1;
    model_num = 1;
    atom_name = "X";
    residue_name = "X";
    chain_name = chain_names[0];
    atom_name_label = "X";
    x = 0;
    y = 0;
    z = 0;
    a = 1;
    b = 0;
  }

  void bind_stream(std::ostream &s) { stream.rdbuf(s.rdbuf()); }

  void write() {
    std::replace(atom_name.begin(), atom_name.end(), '*', '\'');
    if (atom_name == "O1P")
      atom_name = "OP1";
    if (atom_name == "O2P")
      atom_name = "OP2";
    if (std::count_if(chain_name.begin(), chain_name.end(),
                      [](char c) { return c != ' '; }) == 0) {
      chain_name = "A";
    }
    stream << string_format(
                  "ATOM%7d  "
                  "%-4.4s%3.3s%2.2s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%12.12s  ",
                  atom_num, atom_name, residue_name, chain_name, residue_num, x,
                  y, z, a, b, atom_name_label)
           << std::endl;
  }

  template <typename Atom_> void read_atom(const Atom_ &atom) {
    atom_name = atom.name;
    x = atom[0];
    y = atom[1];
    z = atom[2];
    atom_name_label = atom.name.substr(0, 1);
    b = atom.bfactor;
  }

  template <typename F> void write_model_callback(F &&f) {
    write_model_begin();
    f();
    write_model_end();
    model_num++;
    atom_num = 1;
    residue_name = "X";
  }

  template <typename Atom_> void write_atom(const Atom_ &atom) {
    read_atom(atom);
    write();
    atom_num++;
  }

  template <typename Residue_> void write_residue(const Residue_ &residue) {
    residue_name = residue.name;
    if (residue.num != -1) {
      residue_num = residue.num;
    }
    for (auto &&atom : residue) {
      write_atom(atom);
    }
    residue_num++;
  }

  template <typename Chain_> void write_chain(const Chain_ &chain) {
    chain_name = chain.name;
    for (auto &&residue : chain) {
      write_residue(residue);
    }
    write_chain_end();
    residue_num = 1;
    auto it = std::find(chain_names.begin(), chain_names.end(), chain_name);
    chain_name =
        ((it == chain_names.end() || std::next(it) == chain_names.end())
             ? chain_names[0]
             : (*std::next(it)));
  }

  template <typename Model_> void write_model(const Model_ &model) {
    write_model_callback([&]() {
      for (auto &&chain : model) {
        this->write_chain(chain);
      }
    });
  }

  template <typename Pdb_> void write_pdb(const Pdb_ &mol) {
    for (auto &&model : mol) {
      write_model(model);
    }
    write_file_end();
  }

  void write_model_begin() {
    stream << std::left << std::setw(13) << "MODEL" << model_num << std::right
           << std::endl;
  }

  void write_model_end() { stream << "ENDMDL" << std::endl; }

  void write_file_end() { stream << "END" << std::endl; }

  void write_chain_end() {
    stream << "TER " << std::setw(7) << atom_num << "  " << std::left
           << std::setw(4) << " " << std::right << std::setw(3) << residue_name
           << std::setw(2) << chain_name << std::setw(4) << residue_num - 1
           << std::endl;
    atom_num++;
  }
};

inline std::ostream &operator<<(std::ostream &output, const PdbAtom &atom) {
  PdbWriter(output).write_atom(atom);
  return output;
}

inline std::ostream &operator<<(std::ostream &output, const PdbResidue &residue) {
  PdbWriter(output).write_residue(residue);
  return output;
}

inline std::ostream &operator<<(std::ostream &output, const PdbChain &chain) {
  PdbWriter l(output);
  l.write_model_begin();
  l.write_chain(chain);
  l.write_model_end();
  l.write_file_end();
  return output;
}

inline std::ostream &operator<<(std::ostream &output, const PdbModel &model) {
  PdbWriter l(output);
  l.write_model(model);
  l.write_file_end();
  return output;
}

inline std::ostream &operator<<(std::ostream &output, const Pdb &pdb) {
  PdbWriter pdb_writer(output);
  pdb_writer.write_pdb(pdb);
  return output;
}

template<typename Model_>
std::string model_seq(const Model_ &m) {
  auto rs = m.presidues();
  int n = rs.size();
  std::string seq(n, 'X');
  for (int i = 0; i < n; i++) {
    seq[i] = jnc::string_trim_c(rs[i]->name)[0];
  }
  return seq;
}


} // namespace bio
} // namespace jnc


