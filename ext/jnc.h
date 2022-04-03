#pragma once

#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cfloat>

#include <algorithm>
#include <utility>
#include <type_traits>
#include <functional>
#include <memory>

#include <vector>
#include <array>
#include <string>
#include <list>
#include <deque>
#include <array>
#include <set>
#include <tuple>
#include <map>

#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/SVD"
#include "Eigen/Geometry"

namespace jnc {

template<typename T> using Vec = std::vector<T>;
using Str = std::string;

enum {
    JN_CONTINUE,
    JN_BREAK,
    JN_GO
};

template<typename Func_>
void file_each_line(Str fn, Func_ &&f) {
    std::ifstream ifile(fn.c_str());
    Str line;
    while (ifile) {
        std::getline(ifile, line);
        auto r = f(line);
        if (r == JN_CONTINUE)  continue;
        else if (r == JN_BREAK) break;
    }
    ifile.close();
}

template<typename T, typename U>
T lexical_cast(U && u) {
    std::stringstream stream;
    T t;

    stream << u;
    stream >> t;
    return t;
}

#define JN_INT(a) jnc::lexical_cast<int>(a)
#define JN_DBL(a) jnc::lexical_cast<double>(a)
#define JN_FLT(a) jnc::lexical_cast<float>(a)
#define JN_STR(a) jnc::lexical_cast<Str>(a)

template<typename T>
void hash_combine_(std::size_t& seed, const T &value) {
    seed ^= (std::size_t)value + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

template <class T>
inline void hash_combine(std::size_t& seed, const T& v)
{
    std::hash<T> hasher;
    hash_combine_(seed, hasher(v));
}

template <class It>
std::size_t hash_range(It first, It last) {
    std::size_t seed = 0;
    for(; first != last; ++first)
    {
        hash_combine_(seed, *first);
    }
    return seed;
}

template<typename Stream_>
void stream_push(Stream_ && stream) {
}

template<typename Stream_, typename Value_, typename... Values_>
void stream_push(Stream_ && stream, Value_ && value, Values_ && ...values) {
    stream << value;
    stream_push(stream, values...);
}

inline Vec<Str> string_tokenize(const Str &str, const Str &delimiters) {
    Vec<Str> tokens;
    auto lastPos = str.find_first_not_of(delimiters, 0);
    auto pos = str.find_first_of(delimiters, lastPos);
    while (Str::npos != pos || Str::npos != lastPos) {
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiters, pos);
        pos = str.find_first_of(delimiters, lastPos);
    }
    return std::move(tokens);
}

inline Vec<Str> string_tokenize(const Str &str, const Str &delimiters, const Str &temp) {
    Vec<Str> tokens;
    using pair_t = std::pair<Str::size_type, Str::size_type>;
    Vec<pair_t> vec;
    Str::size_type first_i, first_j, second_i, second_j;
    int expected = 0;
    for (Str::size_type i = 0; i < str.size(); i++) {
        int flag = 0;
        Str::size_type j;
        for (j = 0; j < temp.size(); j++) {
            if (str[i] == temp[j]) {
                if (j % 2 == 0 && expected == 0) { flag = 1; break;
                } else if (j % 2 == 1 && expected == 1) { flag = 2; break; }
            }
        }
        if (flag == 1) {
            first_i = i; first_j = j; expected = 1;
        } else if (flag == 2 && j - first_j == 1) {
            second_i = i; second_j = j; expected = 0;
            vec.push_back(std::make_pair(first_i, second_i));
        }
    }
    auto lastPos = str.find_first_not_of(delimiters, 0);
    auto pos = str.find_first_of(delimiters, lastPos);
    while (std::any_of(vec.begin(), vec.end(), [&pos](const pair_t &p){
        return pos != Str::npos && p.first < pos && pos < p.second;
    })) {
        pos = str.find_first_of(delimiters, pos + 1);
    }
    while (Str::npos != pos || Str::npos != lastPos) {
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiters, pos);
        pos = str.find_first_of(delimiters, lastPos);
        while (std::any_of(vec.begin(), vec.end(), [&pos](const pair_t &p){
            return pos != Str::npos && p.first < pos && pos < p.second;
        })) {
            pos = str.find_first_of(delimiters, pos + 1);
        }
    }
    return std::move(tokens);
}

inline void string_upper(Str &s) {
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
}

inline Str string_upper_c(const Str &str) {
    Str s = str; 
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    return s;
}

inline void string_lower(Str &s) {
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
}

inline Str string_lower_c(const Str &str) {
    Str s = str; 
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
}

inline bool string_includes(const Str &s1, const Str &s2) {
    return s1.find(s2) != Str::npos;
}

inline void string_trim(Str &s) {
    if (!(s.empty())) {  
        s.erase(0, s.find_first_not_of(" "));  
        s.erase(s.find_last_not_of(" ") + 1);  
    }
}  

inline Str string_trim_c(const Str &s) {
    if (s.empty()) {
        return s;
    }
    Str::size_type beg = s.find_first_not_of(" \t\n\r");
    Str::size_type end = s.find_last_not_of(" \t\n\r");
    if (beg == Str::npos || end == Str::npos) {
        return "";
    } else {
        return Str(s.begin() + beg, s.begin() + end + 1);
    }
}

inline bool string_starts_with(const Str &str1, const Str &str2) {
    int sz = str2.size();
    for (int i = 0; i < sz; i++) {
        if (str1[i] != str2[i]) return false;
    }
    return true;
}

inline bool string_ends_with(const Str &str1, const Str &str2) {
    int sz = str2.size();
    auto it1 = str1.rbegin();
    auto it2 = str2.rbegin();
    for (int i = 0; i < sz; i++) {
        if (*it1 != *it2) return false;
        it1++;
        it2++;
    }
    return true;
}

template<typename... Strs_>
Str string_merge(Strs_ && ...strs) {
    std::stringstream stream;
    stream_push(stream, strs...);
    return stream.str();
}

///////////////// string_format //////////////////
template<typename T>
T string_to_chars(T &&t) {
    return t;
}

inline const char *string_to_chars(Str &s) {
    return s.c_str();
}

inline const char *string_to_chars(const Str &s) {
    return s.c_str();
}

inline const char *string_to_chars(Str &&s) {
    return s.c_str();
}

template<typename... Pars_>
Str string_format(Str fmt, Pars_ && ...pars) {
    int count = snprintf(NULL, 0, fmt.c_str(), string_to_chars(std::forward<Pars_>(pars))...);
    Str buf;
    buf.resize(count);
    sprintf(&(buf[0]), fmt.c_str(), string_to_chars(pars)...);
    return std::move(buf);
}
//////////////////////////////////////////////////

template<typename _Interval, typename _Ls>
static Str string_join(_Interval && interval, _Ls && ls) {
    std::stringstream stream;
    int i = 0;
    for (const auto & s : ls) {
        if (i != 0) stream << interval;
        stream << s;
        i++;
    }
    return stream.str();
}

#define JN_DIE(...) do {\
    std::cerr << ::jnc::string_merge(__FILE__, "(", __LINE__, "): ", __VA_ARGS__) << std::endl;\
    std::abort();\
} while(0)

#define JN_ASSERT(condition,what) if(!condition){JN_DIE(what);}

namespace pdb {

class Atom : public std::array<double, 3> {
public:
    Str name;
    Str type;
    Str element;
    int num;
    double charge;
    double bfactor = 0;
    bool is_std = true;
};

class Residue : public Vec<Atom> {
public:
    Str name;
    int num = -1;
    bool is_std = true;

    Atom &atom(Str name) {
        for (auto &&atom : *this) {
            if (atom.name == name) return atom;
        }
        JN_DIE("Atom '" + name + "' not found");
    }

    const Atom &atom(Str name) const {
        for (auto &&atom : *this) {
            if (atom.name == name) return atom;
        }
        JN_DIE("Atom '" + name + "' not found");
    }

    bool has_atom(Str name) const {
        for (auto &&atom : *this) {
            if (atom.name == name) return true;
        }
        return false;
    }
};

class Chain : public Vec<Residue> {
public:
    Str name;

    Chain() {}

    Vec<const Atom *> patoms() const {
        Vec<const Atom *> as;
        for (auto && res : *this) {
            for (auto && atom : res) {
                as.push_back(&atom);
            }
        }
        return std::move(as);
    }

    Vec<Atom *> patoms() {
        Vec<Atom *> as;
        for (auto && res : *this) {
            for (auto && atom : res) {
                as.push_back(&atom);
            }
        }
        return std::move(as);
    }
};

class Model : public Vec<Chain> {
public:
    Str name;
    int num;

    Vec<const Atom *> patoms() const {
        Vec<const Atom *> as;
        for (auto && chain : *this) {
            for (auto && res : chain) {
                for (auto && atom : res) {
                    as.push_back(&atom);
                }
            }
        }
        return std::move(as);
    }

    Vec<Atom *> patoms() {
        Vec<Atom *> as;
        for (auto && chain : *this) {
            for (auto && res : chain) {
                for (auto && atom : res) {
                    as.push_back(&atom);
                }
            }
        }
        return std::move(as);
    }

    Vec<const Residue *> presidues() const {
        Vec<const Residue *> rs;
        for (auto && chain : *this) {
            for (auto && res : chain) {
                rs.push_back(&res);
            }
        }
        return std::move(rs);
    }

    Vec<Residue *> presidues() {
        Vec<Residue *> rs;
        int ir = 0;
        for (auto && chain : *this) {
            for (auto && res : chain) {
                if (ir == 0) {
                    std::cout << "presidues" << std::endl;
                    std::cout << &res << std::endl;
                }
                rs.push_back(&res);
                ir++;
            }
        }
        return std::move(rs);
    }
};

class Pdb : public Vec<Model> {
public:
    Str name;

    Pdb() {}

    void remove_hydrogens() {
        for (auto && model : *this) {
            for (auto && chain : model) {
                for (auto && res : chain) {
                    auto r = res;
                    r.clear();
                    for (auto && atom : res) {
                        //                    if (atom.name[0] != 'H' && (!std::isdigit(atom.name[0]) || atom.name[1] != 'H')) 
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
        for (auto && model : *this) {
            for (auto && chain : model) {
                for (auto && res : chain) {
                    std::sort(res.begin(), res.end(), [](const Atom &a1, const Atom &a2){
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
inline void write_line(
    std::ostream &out,
    double x, double y, double z,
    const Str &atom_name, int atom_index,
    const Str &residue_name, int residue_index, const Str &chain_name)
{
    out
        << string_format("ATOM%7d  %-4.4s%3.3s%2.2s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%12.12s  ",
                         atom_index + 1, atom_name, residue_name, chain_name, residue_index + 1, x, y, z, 1.00, 0.00, "" + Str(1, atom_name[0]))
        << std::endl;
}

class PdbReader {
public:
    Residue atoms;
    Chain residues;
    Model chains;
    Pdb &models;

    PdbReader(Pdb &pdb) : models(pdb) {}

    int model_num = 0;

    struct ParsedLine {
        Str atom_name, atom_type, atom_flag, atom_element, res_name, res_flag, chain_name;
        int atom_num, res_num;
        double x, y, z, atom_bfactor;
        bool is_std;
    };

    ParsedLine ol;

    void parse_line(Str line, ParsedLine &pl) {
        pl.res_name     = string_trim_c(line.substr(17, 3));
        pl.res_num      = JN_INT(string_trim_c(line.substr(22, 4)));
        pl.res_flag     = string_trim_c(line.substr(26, 1));
        pl.chain_name   = string_trim_c(line.substr(20, 2));
        pl.is_std       = (!line.compare(0, 4, "ATOM"));

        pl.atom_name    = string_trim_c(line.substr(12, 4));
        pl.atom_num     = JN_INT(string_trim_c(line.substr(6, 5)));
        pl.atom_flag    = string_trim_c(line.substr(16, 1));

        if (line.size() >= 66) {
            pl.atom_bfactor = JN_DBL(line.substr(60, 6)); 
        }
        else {
            pl.atom_bfactor = 0.0;
        }

        if (line.size() >= 78) {
            pl.atom_element = line.substr(77, 1); 
        }
        else {
            pl.atom_element = "X";
        }

        pl.x            = JN_DBL(string_trim_c(line.substr(30, 8)));
        pl.y            = JN_DBL(string_trim_c(line.substr(38, 8)));
        pl.z            = JN_DBL(string_trim_c(line.substr(46, 8)));

        if (line.size() >= 78) {
            pl.atom_type    = string_trim_c(line.substr(76, 2));
        }
        else {
            pl.atom_type    = "X";
        }


        std::replace(pl.atom_name.begin(), pl.atom_name.end(), '*', '\'');
        if (pl.atom_name == "O1P") pl.atom_name = "OP1";
        if (pl.atom_name == "O2P") pl.atom_name = "OP2";
    };

    void add_residue() {
        if (!atoms.empty()) {
            Residue residue = std::move(atoms);
            residue.is_std = std::all_of(atoms.begin(), atoms.end(), [](const Atom &atom){ return atom.is_std; });
            residue.name = ol.res_name;
            residue.num = ol.res_num;
            residues.push_back(std::move(residue));
        }
    };

    void add_chain() {
        add_residue();
        if (!residues.empty()) {
            Chain chain = std::move(residues);
            chain.name = ol.chain_name;
            chains.push_back(std::move(chain));
        }
    };

    void add_atom(Str line) {
        ParsedLine pl;
        parse_line(line, pl);

        if (!atoms.empty()) {
            if (pl.res_num != ol.res_num || pl.res_name != ol.res_name || pl.res_flag != ol.res_flag || pl.chain_name != ol.chain_name) {
                add_residue();
            }
            if (pl.chain_name != ol.chain_name) {
                add_chain();
            }
        }

//        if (string_starts_with(pl.atom_name, "H")) return;
        if (pl.res_name == "HOH" || pl.res_name == "H2O") return;

        // Check whether the atom has appeared before.
        if (!pl.atom_flag.empty()) {
            if (std::any_of(atoms.begin(), atoms.end(), [&pl](const Atom &atom){
                return atom.name == pl.atom_name;
            })) return;
        }

        Atom atom;
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
            Model model = std::move(chains);
            model.num = model_num;
            models.push_back(std::move(model));
            model_num++;
        }
    };

    void read(const Str &fn) {
        std::ifstream ifile(fn.c_str());
        while (ifile) {
            Str line;
            std::getline(ifile, line);
            if (string_starts_with(line, "ATOM")) {
                add_atom(line);
            }
            else if (string_starts_with(line, "HETATM")) {
                add_atom(line);
            }
            else if (string_starts_with(line, "MODEL")) {
                add_model();

                auto &&v = string_tokenize(line, " ");
                if (v.size() >= 2) {
                    model_num = JN_INT(v[1])-1;
                }
            }
            else if (string_starts_with(line, "ENDMDL")) {
                add_model();
            }
            else if (string_starts_with(line, "TER")) {
                add_chain();
            }
            else if (string_trim_c(line) == "END") {
                break;
            }
            else {
                continue;
            }
        }
        ifile.close();

        add_model();

        models.name = fn;
    }
};

Pdb read_pdb(const Str &filename) {
    Pdb pdb;
    PdbReader pdb_reader(pdb);
    pdb_reader.read(filename);
    return std::move(pdb);
}

class PdbWriter {
public:
    const Vec<Str> chain_names = {
        "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N",
        "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"
    };

    int atom_num, residue_num, model_num;
    Str atom_name, residue_name, chain_name, atom_name_label;
    double x, y, z, a, b;

    std::ostream stream;

    PdbWriter() : stream(std::cout.rdbuf()) {
        init();
    }

    PdbWriter(std::ostream &out) : stream(out.rdbuf()) {
        init();
    }

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

    void bind_stream(std::ostream &s) {
        stream.rdbuf(s.rdbuf());
    }

    void write() {
        std::replace(atom_name.begin(), atom_name.end(), '*', '\'');
        if (atom_name == "O1P") atom_name = "OP1";
        if (atom_name == "O2P") atom_name = "OP2";
        if (std::count_if(chain_name.begin(), chain_name.end(), [](char c) {return c != ' '; }) == 0) {
            chain_name = "A";
        }
        stream
            << string_format("ATOM%7d  %-4.4s%3.3s%2.2s%4d%12.3f%8.3f%8.3f%6.2f%6.2f%12.12s  ",
                             atom_num, atom_name, residue_name, chain_name, residue_num, x, y, z, a, b, atom_name_label)
            << std::endl;
    }

    template<typename Atom_>
    void read_atom(const Atom_ &atom) {
        atom_name = atom.name;
        x = atom[0];
        y = atom[1];
        z = atom[2];
        atom_name_label = atom.name.substr(0, 1);
        b = atom.bfactor;
    }

    template<typename F>
    void write_model_callback(F &&f) {
        write_model_begin();
        f();
        write_model_end();
        model_num++;
        atom_num = 1;
        residue_name = "X";
    }

    template<typename Atom_>
    void write_atom(const Atom_ &atom) {
        read_atom(atom);
        write();
        atom_num++;
    }

    template<typename Residue_>
    void write_residue(const Residue_ &residue) {
        residue_name = residue.name;
        if (residue.num != -1) {
            residue_num = residue.num;
        }
        for (auto && atom : residue) {
            write_atom(atom);
        }
        residue_num++;
    }

    template<typename Chain_>
    void write_chain(const Chain_ &chain) {
        chain_name = chain.name;
        for (auto &&residue : chain) {
            write_residue(residue);
        }
        write_chain_end();
        residue_num = 1;
        auto it = std::find(chain_names.begin(), chain_names.end(), chain_name);
        chain_name = ((it == chain_names.end() || std::next(it) == chain_names.end()) ? chain_names[0] : (*std::next(it)));
    }

    template<typename Model_>
    void write_model(const Model_ &model) {
        write_model_callback([&]() {
            for (auto && chain : model) {
                this->write_chain(chain);
            }
        });
    }

    template<typename Pdb_>
    void write_pdb(const Pdb_ &mol) {
        for (auto && model : mol) {
            write_model(model);
        }
        write_file_end();
    }

    void write_model_begin() {
        stream
            << std::left
            << std::setw(13) << "MODEL"
            << model_num
            << std::right
            << std::endl;
    }

    void write_model_end() {
        stream << "ENDMDL" << std::endl;
    }

    void write_file_end() {
        stream << "END" << std::endl;
    }

    void write_chain_end() {
        stream
            << "TER "
            << std::setw(7) << atom_num
            << "  "
            << std::left
            << std::setw(4) << " "
            << std::right
            << std::setw(3) << residue_name
            << std::setw(2) << chain_name
            << std::setw(4) << residue_num - 1
            << std::endl;
        atom_num++;
    }

};

inline std::ostream &operator <<(std::ostream &output, const Atom &atom) {
    PdbWriter(output).write_atom(atom);
    return output;
}

inline std::ostream &operator <<(std::ostream &output, const Residue &residue) {
    PdbWriter(output).write_residue(residue);
    return output;
}

inline std::ostream &operator <<(std::ostream &output, const Chain &chain) {
    PdbWriter l(output);
    l.write_model_begin();
    l.write_chain(chain);
    l.write_model_end();
    l.write_file_end();
    return output;
}

inline std::ostream &operator <<(std::ostream &output, const Model &model) {
    PdbWriter l(output);
    l.write_model(model);
    l.write_file_end();
    return output;
}

inline std::ostream &operator <<(std::ostream &output, const Pdb &pdb) {
    PdbWriter pdb_writer(output);
    pdb_writer.write_pdb(pdb);
    return output;
}

} // namespace pdb

namespace mol2 {

#define MOL2_MOL_TYPE_MAP \
    ELT(SMALL        ,/*MOL_TYPE_*/SMALL       )SEP \
    ELT(BIOPOLYMER   ,/*MOL_TYPE_*/BIOPOLYMER  )SEP \
    ELT(PROTEIN      ,/*MOL_TYPE_*/PROTEIN     )SEP \
    ELT(NUCLEIC_ACID ,/*MOL_TYPE_*/NUCLEIC_ACID)SEP \
    ELT(SACCHARIDE   ,/*MOL_TYPE_*/SACCHARIDE  )SEP 

#define MOL2_CHARGE_TYPE_MAP \
    ELT(NO_CHARGES      ,/*CHARGE_TYPE_*/NO_CHARGES      )SEP \
    ELT(DEL_RE          ,/*CHARGE_TYPE_*/DEL_RE          )SEP \
    ELT(GASTEIGER       ,/*CHARGE_TYPE_*/GASTEIGER       )SEP \
    ELT(GAST_HUCK       ,/*CHARGE_TYPE_*/GAST_HUCK       )SEP \
    ELT(HUCKEL          ,/*CHARGE_TYPE_*/HUCKEL          )SEP \
    ELT(PULLMAN         ,/*CHARGE_TYPE_*/PULLMAN         )SEP \
    ELT(GAUSS80_CHARGES ,/*CHARGE_TYPE_*/GAUSS80_CHARGES )SEP \
    ELT(AMPAC_CHARGES   ,/*CHARGE_TYPE_*/AMPAC_CHARGES   )SEP \
    ELT(MULLIKEN_CHARGES,/*CHARGE_TYPE_*/MULLIKEN_CHARGES)SEP \
    ELT(DICT_CHARGES    ,/*CHARGE_TYPE_*/DICT_CHARGES   )SEP \
    ELT(MMFF94_CHARGES  ,/*CHARGE_TYPE_*/MMFF94_CHARGES  )SEP \
    ELT(USER_CHARGES    ,/*CHARGE_TYPE_*/USER_CHARGES    )SEP 

#define MOL2_BOND_TYPE_MAP \
    ELT("1"  ,/*BOND_TYPE_*/SINGLE       ,single       )SEP \
    ELT("2"  ,/*BOND_TYPE_*/DOUBLE       ,double       )SEP \
    ELT("3"  ,/*BOND_TYPE_*/TRIPLE       ,triple       )SEP \
    ELT("am" ,/*BOND_TYPE_*/AMIDE        ,amide        )SEP \
    ELT("ar" ,/*BOND_TYPE_*/AROMATIC     ,aromatic     )SEP \
    ELT("du" ,/*BOND_TYPE_*/DUMMY        ,dummy        )SEP \
    ELT("un" ,/*BOND_TYPE_*/UNKNOWN      ,unknown      )SEP \
    ELT("nc" ,/*BOND_TYPE_*/NOTCONNECTED ,not connected)SEP

#define MOL2_BOND_STATUS_MAP \
    ELT(1   ,/*BOND_STATUS_*/TYPECOL      )SEP \
    ELT(2   ,/*BOND_STATUS_*/GROUP        )SEP \
    ELT(4   ,/*BOND_STATUS_*/CAP          )SEP \
    ELT(8   ,/*BOND_STATUS_*/BACKBONE     )SEP \
    ELT(16  ,/*BOND_STATUS_*/DICT         )SEP \
    ELT(32  ,/*BOND_STATUS_*/INTERRES     )SEP

#define MOL2_ATOM_TYPE_MAP \
    ELT(Any     ,/*ATOM_TYPE_*/ANY     ,"any atom"                                                                        )SEP \
    ELT(Hev     ,/*ATOM_TYPE_*/HEV     ,"heavy atom (non hydrogen)"                                                       )SEP \
    ELT(Het     ,/*ATOM_TYPE_*/HET     ,"heteroatom = N, O, S, P"                                                         )SEP \
    ELT(LP      ,/*ATOM_TYPE_*/LP      ,"lone pair"                                                                       )SEP \
    ELT(Du      ,/*ATOM_TYPE_*/DU      ,"dummy atom"                                                                      )SEP \
    ELT(Du.C    ,/*ATOM_TYPE_*/DU_C    ,"dummy carbon"                                                                    )SEP \
\
    ELT(H       ,/*ATOM_TYPE_*/H       ,"hydrogen"                                                                        )SEP \
    ELT(H.spc   ,/*ATOM_TYPE_*/H_SPC   ,"hydrogen in Single Point Charge (SPC) water model"                               )SEP \
    ELT(H.t3p   ,/*ATOM_TYPE_*/H_T3P   ,"hydrogen in Transferable intermolecular Potential (TIP3P) water model"           )SEP \
\
    ELT(C.3     ,/*ATOM_TYPE_*/C_3     ,"carbon sp3"                                                                      )SEP \
    ELT(C.2     ,/*ATOM_TYPE_*/C_2     ,"carbon sp2"                                                                      )SEP \
    ELT(C.1     ,/*ATOM_TYPE_*/C_1     ,"carbon sp"                                                                       )SEP \
    ELT(C.ar    ,/*ATOM_TYPE_*/C_AR    ,"carbon aromatic"                                                                 )SEP \
    ELT(C.cat   ,/*ATOM_TYPE_*/C_CAT   ,"carbocation (C +) used only in a guadinium group"                                )SEP \
\
    ELT(N.4     ,/*ATOM_TYPE_*/N_4     ,"nitrogen sp3 positively charged"                                                 )SEP \
    ELT(N.3     ,/*ATOM_TYPE_*/N_3     ,"nitrogen sp3"                                                                    )SEP \
    ELT(N.2     ,/*ATOM_TYPE_*/N_2     ,"nitrogen sp2 "                                                                   )SEP \
    ELT(N.1     ,/*ATOM_TYPE_*/N_1     ,"nitrogen sp"                                                                     )SEP \
    ELT(N.ar    ,/*ATOM_TYPE_*/N_AR    ,"nitrogen aromatic"                                                               )SEP \
    ELT(N.am    ,/*ATOM_TYPE_*/N_AM    ,"nitrogen amide"                                                                  )SEP \
    ELT(N.pl3   ,/*ATOM_TYPE_*/N_PL3   ,"nitrogen trigonal planar"                                                        )SEP \
\
    ELT(O.3     ,/*ATOM_TYPE_*/O_3     ,"oxygen sp3"                                                                      )SEP \
    ELT(O.2     ,/*ATOM_TYPE_*/O_2     ,"oxygen sp2"                                                                      )SEP \
    ELT(O.co2   ,/*ATOM_TYPE_*/O_CO2   ,"oxygen in carboxylate and phosphate groups"                                      )SEP \
    ELT(O.spc   ,/*ATOM_TYPE_*/O_SPC   ,"oxygen in Single Point Charge (SPC) water model "                                )SEP \
    ELT(O.t3p   ,/*ATOM_TYPE_*/O_T3P   ,"oxygen in Transferable Intermolecular Potential (TIP3P) water model"             )SEP \
\
    ELT(S.3     ,/*ATOM_TYPE_*/S_3     ,"sulfur sp3"                                                                      )SEP \
    ELT(S.2     ,/*ATOM_TYPE_*/S_2     ,"sulfur sp2"                                                                      )SEP \
    ELT(S.O     ,/*ATOM_TYPE_*/S_O     ,"sulfoxide sulfur"                                                                )SEP \
    ELT(S.O2    ,/*ATOM_TYPE_*/S_O2    ,"sulfone sulfur"                                                                  )SEP \
    ELT(P.3     ,/*ATOM_TYPE_*/P_3     ,"phosphorous sp3"                                                                 )SEP \
\
    ELT(Hal     ,/*ATOM_TYPE_*/HAL     ,"halogen: fluorine(F), chlorine(Cl), bromine(Br), iodine(I), astatine(At)"        )SEP \
    ELT(F       ,/*ATOM_TYPE_*/F       ,"fluorine"                                                                        )SEP \
    ELT(Cl      ,/*ATOM_TYPE_*/CL      ,"chlorine"                                                                        )SEP \
    ELT(Br      ,/*ATOM_TYPE_*/BR      ,"bromine"                                                                         )SEP \
    ELT(I       ,/*ATOM_TYPE_*/I       ,"iodine"                                                                          )SEP \
\
    ELT(Li      ,/*ATOM_TYPE_*/LI      ,"lithium"                                                                         )SEP \
    ELT(Na      ,/*ATOM_TYPE_*/NA      ,"sodium"                                                                          )SEP \
    ELT(Mg      ,/*ATOM_TYPE_*/MG      ,"magnesium"                                                                       )SEP \
    ELT(Fe      ,/*ATOM_TYPE_*/FE      ,"iron"                                                                            )SEP \
    ELT(Al      ,/*ATOM_TYPE_*/AL      ,"aluminum"                                                                        )SEP \
    ELT(Si      ,/*ATOM_TYPE_*/SI      ,"silicon"                                                                         )SEP \
    ELT(K       ,/*ATOM_TYPE_*/K       ,"potassium"                                                                       )SEP \
    ELT(Ca      ,/*ATOM_TYPE_*/CA      ,"calcium"                                                                         )SEP \
    ELT(Cr.th   ,/*ATOM_TYPE_*/CR_TH   ,"chromium (tetrahedral)"                                                          )SEP \
    ELT(Cr.oh   ,/*ATOM_TYPE_*/CR_OH   ,"chromium (octahedral)"                                                           )SEP \
    ELT(Co.oh   ,/*ATOM_TYPE_*/CO_OH   ,"cobalt (octahedral)"                                                             )SEP \
    ELT(Mn      ,/*ATOM_TYPE_*/MN      ,"manganese"                                                                       )SEP \
    ELT(Cu      ,/*ATOM_TYPE_*/CU      ,"copper"                                                                          )SEP \
    ELT(Zn      ,/*ATOM_TYPE_*/ZN      ,"zinc"                                                                            )SEP \
    ELT(Se      ,/*ATOM_TYPE_*/SE      ,"selenium"                                                                        )SEP \
    ELT(Mo      ,/*ATOM_TYPE_*/MO      ,"molybdenum"                                                                      )SEP \
    ELT(Sn      ,/*ATOM_TYPE_*/SN      ,"tin"                                                                             )SEP

#define MOL2_ATOM_STATUS_MAP \
    ELT(1   ,/*ATOM_STATUS_*/DSPMOD       )SEP \
    ELT(2   ,/*ATOM_STATUS_*/TYPECOL      )SEP \
    ELT(4   ,/*ATOM_STATUS_*/CAP          )SEP \
    ELT(8   ,/*ATOM_STATUS_*/BACKBONE     )SEP \
    ELT(16  ,/*ATOM_STATUS_*/DICT         )SEP \
    ELT(32  ,/*ATOM_STATUS_*/ESSENTIAL    )SEP \
    ELT(64  ,/*ATOM_STATUS_*/WATER        )SEP \
    ELT(128 ,/*ATOM_STATUS_*/DIRECT       )SEP

#define JN_MOL2_ATOM_TABLE \
    ELT("H"       , ATOM_TYPE_H    , ATOM_ELEMENT_H) SEP \
    ELT("H.SPC"   , ATOM_TYPE_H_SPC, ATOM_ELEMENT_H) SEP \
    ELT("H.T3P"   , ATOM_TYPE_H_T3P, ATOM_ELEMENT_H) SEP \
\
    ELT("C.1"     , ATOM_TYPE_C_1  , ATOM_ELEMENT_C) SEP \
    ELT("C.2"     , ATOM_TYPE_C_2  , ATOM_ELEMENT_C) SEP \
    ELT("C.3"     , ATOM_TYPE_C_3  , ATOM_ELEMENT_C) SEP \
    ELT("C.AR"    , ATOM_TYPE_C_AR , ATOM_ELEMENT_C) SEP \
    ELT("C.CAT"   , ATOM_TYPE_C_CAT, ATOM_ELEMENT_C) SEP \
\
    ELT("N.1"     , ATOM_TYPE_N_1  , ATOM_ELEMENT_N) SEP \
    ELT("N.2"     , ATOM_TYPE_N_2  , ATOM_ELEMENT_N) SEP \
    ELT("N.3"     , ATOM_TYPE_N_3  , ATOM_ELEMENT_N) SEP \
    ELT("N.4"     , ATOM_TYPE_N_4  , ATOM_ELEMENT_N) SEP \
    ELT("N.AR"    , ATOM_TYPE_N_AR , ATOM_ELEMENT_N) SEP \
    ELT("N.AM"    , ATOM_TYPE_N_AM , ATOM_ELEMENT_N) SEP \
    ELT("N.PL3"   , ATOM_TYPE_N_PL3, ATOM_ELEMENT_N) SEP \
\
    ELT("O.2"     , ATOM_TYPE_O_2  , ATOM_ELEMENT_O) SEP \
    ELT("O.3"     , ATOM_TYPE_O_3  , ATOM_ELEMENT_O) SEP \
    ELT("O.CO2"   , ATOM_TYPE_O_CO2, ATOM_ELEMENT_O) SEP \
    ELT("O.T3P"   , ATOM_TYPE_O_T3P, ATOM_ELEMENT_O) SEP \
    ELT("O.SPC"   , ATOM_TYPE_O_SPC, ATOM_ELEMENT_O) SEP \
\
    ELT("P.3"     , ATOM_TYPE_P_3  , ATOM_ELEMENT_P) SEP \
\
    ELT("S.2"     , ATOM_TYPE_S_2  , ATOM_ELEMENT_S) SEP \
    ELT("S.3"     , ATOM_TYPE_S_3  , ATOM_ELEMENT_S) SEP \
    ELT("S.O"     , ATOM_TYPE_S_O  , ATOM_ELEMENT_S) SEP \
    ELT("S.O2"    , ATOM_TYPE_S_O2 , ATOM_ELEMENT_S) SEP \
\
    ELT("K"       , ATOM.TYPE_K    , ATOM_ELEMENT_K ) SEP \
    ELT("F"       , ATOM.TYPE_F    , ATOM_ELEMENT_F ) SEP \
    ELT("I"       , ATOM.TYPE_I    , ATOM_ELEMENT_I ) SEP \
    ELT("LP"      , ATOM.TYPE_LP   , ATOM_ELEMENT_LP) SEP \
    ELT("MG"      , ATOM.TYPE_MG   , ATOM_ELEMENT_MG) SEP \
    ELT("CR.OH"   , ATOM_TYPE_CR_OH, ATOM_ELEMENT_CR) SEP \
    ELT("CR.TH"   , ATOM_TYPE_CR_TH, ATOM_ELEMENT_CR) SEP \
    ELT("CO.OH"   , ATOM_TYPE_CO_OH, ATOM_ELEMENT_CO) SEP \
    ELT("CL"      , ATOM.TYPE_CL   , ATOM_ELEMENT_CL) SEP \
    ELT("SE"      , ATOM.TYPE_SE   , ATOM_ELEMENT_SE) SEP \
    ELT("NA"      , ATOM.TYPE_NA   , ATOM_ELEMENT_NA) SEP \
    ELT("SI"      , ATOM.TYPE_SI   , ATOM_ELEMENT_SI) SEP \
    ELT("FE"      , ATOM.TYPE_FE   , ATOM_ELEMENT_FE) SEP \
    ELT("ZN"      , ATOM.TYPE_ZN   , ATOM_ELEMENT_ZN) SEP \
    ELT("SN"      , ATOM.TYPE_SN   , ATOM_ELEMENT_SN) SEP \
    ELT("LI"      , ATOM.TYPE_LI   , ATOM_ELEMENT_LI) SEP \
    ELT("AL"      , ATOM.TYPE_AL   , ATOM_ELEMENT_AL) SEP \
    ELT("CA"      , ATOM.TYPE_CA   , ATOM_ELEMENT_CA) SEP \
    ELT("MN"      , ATOM.TYPE_MN   , ATOM_ELEMENT_MN) SEP \
    ELT("CU"      , ATOM.TYPE_CU   , ATOM_ELEMENT_CU) SEP \
    ELT("BR"      , ATOM.TYPE_BR   , ATOM_ELEMENT_BR) SEP \
    ELT("MO"      , ATOM.TYPE_MO   , ATOM_ELEMENT_MO) SEP \
\
    ELT("DU"      , ATOM.TYPE_DU   , ATOM_ELEMENT_X) SEP \
    ELT("ANY"     , ATOM.TYPE_ANY  , ATOM_ELEMENT_X) SEP \
    ELT("DU.C"    , ATOM_TYPE_DU_C , ATOM_ELEMENT_X) SEP \
    ELT("HAL"     , ATOM.TYPE_HAL  , ATOM_ELEMENT_X) SEP \
    ELT("HET"     , ATOM.TYPE_HET  , ATOM_ELEMENT_X) SEP \
    ELT("HEV"     , ATOM.TYPE_HEV  , ATOM_ELEMENT_X)

/**
 * MolType
 */
class MolType {
public:
    enum {
#define ELT(name, type) type
#define SEP ,
        MOL2_MOL_TYPE_MAP
#undef ELT
#undef SEP
    };

    static Str get_name(int id) {
        std::map<int, Str> m {
#define ELT(name, type) {type, #name}
#define SEP ,
            MOL2_MOL_TYPE_MAP
#undef ELT
#undef SEP
        };
        return m[id];
    }

    static int get_id(const Str &name) {
        static std::map<Str, int> m {
#define ELT(name, type) {#name, type}
#define SEP ,
            MOL2_MOL_TYPE_MAP
#undef ELT
#undef SEP
        };
        return m[name];
    }
};

/**
 * Charge Type
 */
class ChargeType {
public:
    enum {
#define ELT(name, type) type
#define SEP ,
        MOL2_CHARGE_TYPE_MAP
#undef ELT
#undef SEP
    };

    static Str get_name(int id) {
        std::map<int, Str> m {
#define ELT(name, type) {type, #name}
#define SEP ,
            MOL2_CHARGE_TYPE_MAP
#undef ELT
#undef SEP
        };
        return m[id];
    }

    static int get_id(const Str &name) {
        std::map<Str, int> m {
#define ELT(name, type) {#name, type}
#define SEP ,
            MOL2_CHARGE_TYPE_MAP
#undef ELT
#undef SEP
        };
        return m[name];
    }
};

/**
 * Bond Type
 */
class BondType {
public:
    enum {
#define ELT(name, type, annotation) type
#define SEP ,
        MOL2_BOND_TYPE_MAP
#undef ELT
#undef SEP
    };

    static Str get_name(int id) {
        std::map<int, Str> m {
#define ELT(name, type, annotation) {type, name}
#define SEP ,
            MOL2_BOND_TYPE_MAP
#undef ELT
#undef SEP
        };
        return m[id];
    }

    static int get_id(const Str &name) {
        std::map<Str, int> m {
#define ELT(name, type, annotation) {name, type}
#define SEP ,
            MOL2_BOND_TYPE_MAP
#undef ELT
#undef SEP
        };
        return m[name];
    }
};

/**
 * BondStatus
 */
class BondStatus {
public:
    enum {
#define ELT(n, status) status = n
#define SEP ,
        MOL2_BOND_STATUS_MAP
#undef ELT
#undef SEP
    };

    static Str get_name(int id) {
        std::map<int, Str> m {
#define ELT(n, status) {status, #status}
#define SEP ,
            MOL2_BOND_STATUS_MAP
#undef ELT
#undef SEP
        };
        return m[id];

    }

    static int get_id(const Str &name) {
        std::map<Str, int> m {
#define ELT(n, status) {#status, status}
#define SEP ,
            MOL2_BOND_STATUS_MAP
#undef ELT
#undef SEP
        };
        return m[name];

    }
};

/**
 * AtomType
 */
class AtomType {
public:
    enum {
#define ELT(name, type, annotation) type
#define SEP ,
        MOL2_ATOM_TYPE_MAP
#undef ELT
#undef SEP
    };

    static Str get_name(int id) {
        std::map<int, Str> m {
#define ELT(name, type, annotation) {type, #name}
#define SEP ,
            MOL2_ATOM_TYPE_MAP
#undef ELT
#undef SEP
        };
        return m[id];

    }

    static int get_id(const Str &name) {
        std::map<Str, int> m {
#define ELT(name, type, annotation) {#name, type}
#define SEP ,
            MOL2_ATOM_TYPE_MAP
#undef ELT
#undef SEP
        };
        return m[name];
    }
};

/**
 * AtomStatus related definitions
 */
class AtomStatus {
public:
    enum {
#define ELT(n, status) status = n
#define SEP ,
        MOL2_ATOM_STATUS_MAP
#undef ELT
#undef SEP
    };

    static Str get_name(int id) {
        std::map<int, Str> m {
#define ELT(n, status) {status, #status}
#define SEP ,
            MOL2_ATOM_STATUS_MAP
#undef ELT
#undef SEP
        };
        return m[id];

    }

    static int get_id(const Str &name) {
        std::map<Str, int> m {
#define ELT(n, status) {#status, status}
#define SEP ,
            MOL2_ATOM_STATUS_MAP
#undef ELT
#undef SEP
        };
        return m[name];

    }
};

struct Bond {
    int origin_atom_id;
    int target_atom_id;
    Str type;
    Str status_bits = "";
};
using Bonds = Vec<Bond>;

struct Atom {
    Str name;
    double x;
    double y;
    double z;
    Str type;
    int subst_id = INT_MAX;
    Str subst_name = "";
    double charge = DBL_MAX;
    Str status_bit = "";

    double &operator[](int i) {
        return i == 0 ? x : (i == 1 ? y : z);
    }

    const double &operator[](int i) const {
        return i == 0 ? x : (i == 1 ? y : z);
    }
};
using Atoms = Vec<Atom>;

struct Substructure {
    Str name;
    int root_atom;
    Str type;
    int dict_type = INT_MAX;
    Str chain;
    Str sub_type;
    int inter_bonds = INT_MAX;
    Str status;
    Str comment;
};
using Substructures = Vec<Substructure>;

struct Mol2 {
    Str mol_name;

    int num_atoms;
    int num_bonds = INT_MAX;
    int num_subst = INT_MAX;
    int num_feat = INT_MAX;
    int num_sets = INT_MAX;

    Str mol_type;
    Str charge_type;
    Str status_bits = "";
    Str mol_comment = "";

    Atoms atoms;
    Bonds bonds;
    Substructures substructures;

    void remove_hydrogens() {
        Atoms as;
        Vec<int> vi(atoms.size(), -1);
        int i = 0;
        int j = 0;
        for (auto &&atom : atoms) {
            auto &type = atom.type;
            if (type != "H" && type.find("H.") == Str::npos) {
                vi[i] = j;
                as.push_back(atom);
                j++;
            }
            i++;
        }
        atoms = std::move(as);

        Bonds bs;
        for (auto &&bond : bonds) {
            bond.origin_atom_id = vi[bond.origin_atom_id];
            bond.target_atom_id = vi[bond.target_atom_id];
            if (bond.origin_atom_id != -1 && bond.target_atom_id != -1) bs.push_back(bond);
        }
        bonds = std::move(bs);

        num_atoms = atoms.size();
        num_bonds = bonds.size();
    }

    void write(Str fn) {
        std::ofstream ofile(fn.c_str());
        write(ofile);
        ofile.close();
    }

    void write(std::ostream &output) {
        write_molecule(output);
        write_atoms(output);
        write_bonds(output);
        write_substructures(output);
    }

    void write_molecule(std::ostream &out) {
        out << "@<TRIPOS>MOLECULE" << std::endl;
        out << mol_name << std::endl;
        out << num_atoms;
        if (num_bonds != INT_MAX) {
            out << ' ' << num_bonds;
            if (num_subst != INT_MAX) {
                out << ' ' << num_subst;
                if (num_feat != INT_MAX) {
                    out << ' ' << num_feat;
                    if (num_sets != INT_MAX) {
                        out << ' ' << num_sets;
                    }
                }
            }
        }
        out << std::endl;

        out << mol_type << std::endl;
        out << charge_type << std::endl;

        if (!status_bits.empty()) out << status_bits << std::endl;
        if (!mol_comment.empty()) out << mol_comment << std::endl;
    }

    void write_atoms(std::ostream &out) {
        out << "@<TRIPOS>ATOM" << std::endl;
        int i = 0;
        for (auto &&atom : atoms) {
            out << i+1 << ' ' << atom.name << ' ' << atom.x << ' ' << atom.y << ' ' << atom.z << ' ' << atom.type;
            if (atom.subst_id != INT_MAX) {
                out << ' ' << atom.subst_id+1;
                if (!atom.subst_name.empty()) {
                    out << ' ' << atom.subst_name;
                    if (atom.charge != DBL_MAX) {
                        out << ' ' << atom.charge;
                        if (!atom.status_bit.empty()) {
                            out << ' ' << atom.status_bit;
                        }
                    }
                }
            }
            out << std::endl;
            i++;
        }
    }

    void write_bonds(std::ostream &out) {
        if (bonds.empty()) return;

        out << "@<TRIPOS>BOND" << std::endl;
        int i = 0;
        for (auto &&bond : bonds) {
            out << i+1 << ' ' << bond.origin_atom_id+1 << ' ' << bond.target_atom_id+1 << ' ' << bond.type;
            if (!bond.status_bits.empty()) out << ' ' << bond.status_bits;
            out << std::endl;
            i++;
        }
    }

    void write_substructures(std::ostream &out) {
        if (substructures.empty()) return;

        out << "@<TRIPOS>SUBSTRUCTURE" << std::endl;
        int i = 0;
        for (auto &&subst : substructures) {
            out << i+1 << ' ' << subst.name << ' ' << subst.root_atom;
            if (!subst.type.empty()) {
                out << ' ' << subst.type;
                if (subst.dict_type != INT_MAX) {
                    out << ' ' << subst.dict_type;
                    if (!subst.chain.empty()) {
                        out << ' ' << subst.chain;
                        if (!subst.sub_type.empty()) {
                            out << ' ' << subst.sub_type;
                            if (subst.inter_bonds != INT_MAX) {
                                out << ' ' << subst.inter_bonds;
                                if (!subst.status.empty()) {
                                    out << ' ' << subst.status;
                                    if (!subst.comment.empty()) {
                                        out << ' ' << subst.comment;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            out << std::endl;
            i++;
        }
    }

};

Vec<Mol2> read_mol2s(const Str &fn) {
    Str section;
    Vec<Mol2> ms;
    Mol2 *p = NULL;

    int mol_lines = 0;
    file_each_line(fn, [&section, &mol_lines, &p, &ms](Str &line){
        if (!line.compare(0, 9, "@<TRIPOS>")) {
            section = line.substr(9);
            if (section == "MOLECULE") {
                ms.push_back(Mol2{});
                p = &(ms.back());
            }
        }
        else if (string_tokenize(line, " \t").empty()) {
        }
        else {
            if (section == "MOLECULE") {
                mol_lines++;
                if (mol_lines == 1) p->mol_name = string_trim_c(line);
                else if (mol_lines == 2) {
                    auto &&vs = string_tokenize(line, " \t");
                    if (vs.size() >= 1) p->num_atoms = JN_INT(vs[0]);
                    if (vs.size() >= 2) p->num_bonds = JN_INT(vs[1]);
                    if (vs.size() >= 3) p->num_subst = JN_INT(vs[2]);
                    if (vs.size() >= 4) p->num_feat  = JN_INT(vs[3]);
                    if (vs.size() >= 5) p->num_sets  = JN_INT(vs[4]);
                }
                else if (mol_lines == 3) {
                    p->mol_type = string_trim_c(line);
                }
                else if (mol_lines == 4) p->charge_type = string_trim_c(line);
                else if (mol_lines == 5) p->status_bits = string_trim_c(line);
                else if (mol_lines == 6) p->mol_comment = string_trim_c(line);
            }
            else if (section == "ATOM") {
                auto &&vs = string_tokenize(line, " \t");
                Atom atom;
                if (vs.size() >= 6) {
                    atom.name = vs[1];
                    atom.x = JN_DBL(vs[2]);
                    atom.y = JN_DBL(vs[3]);
                    atom.z = JN_DBL(vs[4]);
                    atom.type = vs[5];
                    if (vs.size() >= 7) atom.subst_id = JN_INT(vs[6])-1;
                    if (vs.size() >= 8) atom.subst_name = vs[7];
                    if (vs.size() >= 9) atom.charge = JN_DBL(vs[8]);
                    if (vs.size() >= 10) atom.status_bit = vs[9];

                    p->atoms.push_back(atom);
                }
                else {
                    JN_DIE("Wrong mol2 format!");
                }
            }
            else if (section == "BOND") {
                auto &&vs = string_tokenize(line, " \t");
                if (vs.size() >= 4) {
                    Bond bond;
                    bond.origin_atom_id = JN_INT(vs[1])-1;
                    bond.target_atom_id = JN_INT(vs[2])-1;
                    bond.type = vs[3];
                    if (vs.size() > 4) {
                        bond.status_bits  = vs[4];
                    }
                    p->bonds.push_back(bond);
                }
            }
            else if (section == "SUBSTRUCTURE") {
                auto &&vs = string_tokenize(line, " \t");
                if (vs.size() >= 3) {
                    Substructure subst;
                    subst.name = vs[1];
                    subst.root_atom = JN_INT(vs[2]);
                    if (vs.size() > 3) subst.type = vs[3];
                    if (vs.size() > 4) subst.dict_type = JN_INT(vs[4]);
                    if (vs.size() > 5) subst.chain = vs[5];
                    if (vs.size() > 6) subst.sub_type = vs[6];
                    if (vs.size() > 7) subst.inter_bonds = JN_INT(vs[7]);
                    if (vs.size() > 8) subst.status = vs[8];
                    if (vs.size() > 9) subst.comment = vs[9];

                    p->substructures.push_back(subst);
                }
            }
        }
        return JN_GO;
    });
    return std::move(ms);
}

} // namespace mol2

namespace geom {

template<typename NumType>
using MatX = Eigen::Matrix<NumType, -1, -1>;
using Matd = MatX<double>;

template<typename NumType>
using VecX = Eigen::Matrix<NumType, -1, 1>;
using Vecd = VecX<double>;

template<typename T, typename U>
void translate(T &&t, const U &u) {
    for (int i = 0; i < 3; i++) t[i] += u[i];
}

/**
 * Rotate fixed in the origin
 */
template<typename T, typename _Mat>
void rotate(T &&t, const _Mat &mat) {
    auto x = t[0] * mat(0, 0) + t[1] * mat(1, 0) + t[2] * mat(2, 0);
    auto y = t[0] * mat(0, 1) + t[1] * mat(1, 1) + t[2] * mat(2, 1);
    auto z = t[0] * mat(0, 2) + t[1] * mat(1, 2) + t[2] * mat(2, 2);
    t[0] = x; t[1] = y; t[2] = z;
}

/**
 * Rotate fixed in an specified point
 */
template<typename T, typename U>
void rotate(T &&t, const U &origin, const MatX<double> &mat) {
    for (int i = 0; i < 3; i++) t[i] -= origin[i];
    rotate(t, mat);
    for (int i = 0; i < 3; i++) t[i] += origin[i];
}

template<class C1, class C2>
MatX<double> x_rot_mat(C1 c, C2 s) {
    MatX<double> rot_mat(3, 3);
    rot_mat <<
        1, 0, 0,
        0, c, s,
        0, -s, c;
    return rot_mat;
}

template<class C1, class C2>
MatX<double> y_rot_mat(C1 c, C2 s) {
    MatX<double> rot_mat(3, 3);
    rot_mat <<
        c, 0, -s,
        0, 1, 0,
        s, 0, c;
    return rot_mat;
}

template<class C1, class C2>
MatX<double> z_rot_mat(C1 c, C2 s) {
    MatX<double> rot_mat(3, 3);
    rot_mat <<
        c, s, 0,
        -s, c, 0,
        0, 0, 1;
    return rot_mat;
}

template<typename T>
MatX<double> rot_mat(int i, T &&v) {
    double c = std::cos(v);
    double s = std::sin(v);
    assert(i >= 0 && i < 3);
    return (i == 0 ? x_rot_mat(c, s) : (i == 1 ? y_rot_mat(c, s) : z_rot_mat(c, s)));
}

// Rotate along with an axis
class RotateAlong {
public:
    using mat_t = MatX<double>;
    using vec_t = VecX<double>;

    mat_t m_rm;
    vec_t m_beg, m_end;

    RotateAlong() = default;

    template<typename Begin, typename End>
    RotateAlong(const Begin &begin, const End &end, double angle) {
        init(begin, end, angle);
    }

    template<typename Begin, typename End>
    void init(const Begin &begin, const End &end, double t) {
        double r1, r2, c1, c2, s1, s2;
        //L l = end - begin;
        vec_t l(3);
        for (int i = 0; i < 3; i++) l[i] = end[i] - begin[i];

        m_beg.resize(3);
        m_end.resize(3);
        for (int i = 0; i < 3; i++) {
            m_beg[i] = begin[i];
            m_end[i] = end[i];
        }

        r1 = std::sqrt(l[0] * l[0] + l[1] * l[1]);
        c1 = l[1] / r1;
        s1 = l[0] / r1;

        r2 = std::sqrt(l[0] * l[0] + l[1] * l[1] + l[2] * l[2]);
        c2 = l[2] / r2;
        s2 = std::sqrt(l[0] * l[0] + l[1] * l[1]) / r2;

        m_rm = mat_t::Identity(3, 3);
        if (r1 != 0) m_rm = m_rm * z_rot_mat(c1, s1);
        if (r2 != 0) m_rm = m_rm * x_rot_mat(c2, s2);
        m_rm = m_rm * z_rot_mat(std::cos(t), std::sin(t));
        if (r2 != 0) m_rm = m_rm * x_rot_mat(c2, -s2);
        if (r1 != 0) m_rm = m_rm * z_rot_mat(c1, -s1);

    }

    template<typename T>
    RotateAlong &operator ()(T &&t) {
        rotate(t, m_beg, m_rm);
        return *this;
    }
};

template<typename T>
T square(T n) {
    return n*n;
}

template<class P1, class P2> 
double distance(const P1 &p1, const P2 &p2) {
    return std::sqrt(square(p1[0]-p2[0])+square(p1[1]-p2[1])+square(p1[2]-p2[2]));
}

template<class P1, class P2> 
double dist2(const P1 &p1, const P2 &p2) {
    return square(p1[0]-p2[0])+square(p1[1]-p2[1])+square(p1[2]-p2[2]);
}

template<typename P1, typename P2>
inline double norm(P1 &&p1, P2 &&p2, int n) {
    double sum = 0; for (int i = 0; i < n; i++) sum += (p1[i] - p2[i]) * (p1[i] - p2[i]);
    return std::sqrt(sum);
}

template<class P1, class P2, class P3>
inline double angle(const P1 &p1, const P2 &p2, const P3 &p3) {
    double a1 = p1[0] - p2[0], a2 = p1[1] - p2[1], a3 = p1[2] - p2[2];
    double b1 = p3[0] - p2[0], b2 = p3[1] - p2[1], b3 = p3[2] - p2[2];
    double ra = sqrt(a1 * a1 + a2 * a2 + a3 * a3);
    double rb = sqrt(b1 * b1 + b2 * b2 + b3 * b3);
    return acos((a1 * b1 + a2 * b2 + a3 * b3) / (ra * rb));
}

template<class T1, class T2, class T3, class T4>
inline double chirality(const T1 &p1, const T2 &p2, const T3 &p3, const T4 &p4) {
    double a[3] = {p1[0] - p4[0], p1[1] - p4[1], p1[2] - p4[2]};
    double b[3] = {p2[0] - p4[0], p2[1] - p4[1], p2[2] - p4[2]};
    double c[3] = {p3[0] - p4[0], p3[1] - p4[1], p3[2] - p4[2]};
    double d[3] = {b[1]*c[2]-b[2]*c[1], b[2]*c[0]-b[0]*c[2], b[0]*c[1]-b[1]*c[0]};
    return a[0]*d[0]+a[1]*d[1]+a[2]*d[2];
}

template<class P = Eigen::Vector3d, class P1, class P2, class P3>
P normal_vector(const P1 &p1, const P2 &p2, const P3 &p3) {
    double a1, a2, a3, b1, b2, b3, r;
    P p;
    a1 = p2[0] - p1[0]; a2 = p2[1] - p1[1]; a3 = p2[2] - p1[2];
    b1 = p3[0] - p2[0]; b2 = p3[1] - p2[1]; b3 = p3[2] - p2[2];
    p[0] = a2 * b3 - a3 * b2;
    p[1] = b1 * a3 - a1 * b3;
    p[2] = a1 * b2 - a2 * b1;
    r = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    p[0] /= r; p[1] /= r; p[2] /= r;
    return p;
}

template<class P1, class P2, class P3, class P4>
inline double dihedral(const P1 &p1, const P2 &p2, const P3 &p3, const P4 &p4) {
    double a1, a2, a3, b1, b2, b3, c1, c2, c3;
    double a1_, a2_, a3_, b1_, b2_, b3_, c1_, c2_, c3_;
    double r, c, s;

    a1 = p1[0] - p2[0]; a2 = p1[1] - p2[1]; a3 = p1[2] - p2[2];
    b1 = p3[0] - p2[0]; b2 = p3[1] - p2[1]; b3 = p3[2] - p2[2];
    c1 = p4[0] - p2[0]; c2 = p4[1] - p2[1]; c3 = p4[2] - p2[2];
    if (b1 * b1 + b2 * b2 != 0) {
        r = sqrt(b1 * b1 + b2 * b2);
        c = b1 / r; s = -b2 / r;
        a1_ = c * a1 - s * a2; a2_ = s * a1 + c * a2; a3_ = a3;
        b1_ = c * b1 - s * b2; b2_ = s * b1 + c * b2; b3_ = b3;
        c1_ = c * c1 - s * c2; c2_ = s * c1 + c * c2; c3_ = c3;
    } else {
        a1_ = a1; a2_ = a2; a3_ = a3;
        b1_ = b1; b2_ = b2; b3_ = b3;
        c1_ = c1; c2_ = c2; c3_ = c3;
    }
    if (b1_ * b1_ + b3_ * b3_ != 0) {
        r = sqrt(b1_ * b1_ + b3_ * b3_);
        c = b1_ / r; s = b3_ / r;
        a1 = c * a1_ + s * a3_; a2 = a2_; a3 = -s * a1_ + c * a3_;
        b1 = c * b1_ + s * b3_; b2 = b2_; b3 = -s * b1_ + c * b3_;
        c1 = c * c1_ + s * c3_; c2 = c2_; c3 = -s * c1_ + c * c3_;
    } else {
        a1 = a1_; a2 = a2_; a3 = a3_;
        b1 = b1_; b2 = b2_; b3 = b3_;
        c1 = c1_; c2 = c2_; c3 = c3_;
    }
    if (a2 * a2 + a3 * a3 != 0) {
        r = sqrt(a2 * a2 + a3 * a3);
        c = a3 / r; s = a2 / r;
        a1_ = a1; a2_ = c * a2 - s * a3; a3_ = s * a2 + c * a3;
        b1_ = b1; b2_ = c * b2 - s * b3; b3_ = s * b2 + c * b3;
        c1_ = c1; c2_ = c * c2 - s * c3; c3_ = s * c2 + c * c3;
    }
    if (c2_ * c2_ + c3_ * c3_ == 0) return 0;
    double temp = acos(c3_ / sqrt(c2_ * c2_ + c3_ * c3_));
    if (c2_ > 0) {
        temp = -temp;
    }
    return temp;
} 

template<typename NumType>
class SupPos {
public:
    using mat_t = MatX<NumType>;
    using vec_t = VecX<NumType>;

    mat_t rot;
    vec_t c1;
    vec_t c2;
    NumType rmsd;

    SupPos() = default;

    SupPos(const mat_t &m, const mat_t &n) {
        init(m, n);
    }

    void init(const mat_t &m, const mat_t &n) {
        int i, j, len;
        mat_t x, y, g, u, v, I, d;
        std::ostringstream stream;

        if (m.rows() != n.rows() || m.cols() != 3 || n.cols() != 3) {
            stream << "jian::geom::suppos error! (" << m.rows() << ' ' << m.cols() << ") (" << n.rows() << ' ' << n.cols() << ")\n";
            throw stream.str();
        }

        len = m.rows();
        x = m;
        y = n;
        c1 = vec_t::Zero(3);
        c2 = vec_t::Zero(3);
        for (i = 0; i < len; i++) for (j = 0; j < 3; j++) { c1[j] += x(i, j); c2[j] += y(i, j); }
        for (i = 0; i < 3; i++) { c1[i] = c1[i] / len; c2[i] = c2[i] / len; }
        for (i = 0; i < len; i++) for (j = 0; j < 3; j++) { x(i, j) -= c1[j]; y(i, j) -= c2[j]; }

        g = x.transpose() * y;
        Eigen::JacobiSVD<mat_t> svd(g, Eigen::ComputeFullU | Eigen::ComputeFullV);
        u = svd.matrixU();
        v = svd.matrixV();

        I = mat_t::Identity(3, 3);
        if (g.determinant() < 0) I(2, 2) = -1;
        rot = u * I * v.transpose();

        d = x * rot - y;
        rmsd = 0;
        for (int i = 0; i < len; i++) for (int j = 0; j < 3; j++) rmsd += d(i, j) * d(i, j);
        rmsd = std::sqrt(rmsd / len);

        c1 = -c1;
    }

    template<typename T>
    void apply(T &&t) {
        translate(t, c1);
        rotate(t, rot);
        translate(t, c2);
    }

    template<typename M>
    void apply_m(M &m) {
        int i, j;
        for (i = 0; i < m.rows(); i++) {
            for (j = 0; j < 3; j++) {
                m(i, j) += c1[j];
            }
        }
        m = m * rot;
        for (i = 0; i < m.rows(); i++) {
            for (j = 0; j < 3; j++) {
                m(i, j) += c2[j];
            }
        }
    }
};

template<typename NumType>
SupPos<NumType> suppos(const MatX<NumType> &m, const MatX<NumType> &n) {
    return SupPos<NumType>(m, n);
}

template<typename T, typename NumType>
SupPos<NumType> suppos(T &t, const MatX<NumType> &m, const MatX<NumType> &n) {
    SupPos<NumType> sp(m, n);
    sp.apply_m(t);
    return sp;
}

template<typename NumType, typename U, typename F, typename V>
SupPos<NumType> suppos(MatX<NumType> &src, const U &src_indices, const F &tgt, const V &tgt_indices) {
    MatX<NumType> m(src_indices.size(), 3), n(tgt_indices.size(), 3);
    for (int i = 0; i < src_indices.size(); i++) for (int j = 0; j < 3; j++) {
        m(i, j) = src(src_indices[i], j);
        n(i, j) = tgt(tgt_indices[i], j);
    }

    auto sp = suppos(m, n);

    sp.apply_m(src);
    return sp;
}

template<typename NumType1, typename NumType2>
NumType1 rmsd(const MatX<NumType1> &x, const MatX<NumType2> &y) {
    return SupPos<NumType1>(x, y).rmsd;
}

template<typename NumType, typename A, typename B, typename C, typename D>
MatX<NumType> suppos_axis_polar(A theta_o, B phi_o, C theta_n, D phi_n) {
    MatX<NumType> m = MatX<NumType>::Identity(3, 3);
    double ang = - phi_o;
    if (ang != 0) m *= z_rot_mat(std::cos(ang), std::sin(ang));
    ang = -theta_o + theta_n;
    if (ang != 0) m *= y_rot_mat(std::cos(ang), std::sin(ang));
    ang = phi_n;
    if (ang != 0) m *= z_rot_mat(std::cos(ang), std::sin(ang));
    return m;
}


template<typename NumType, typename O, typename N>
MatX<NumType> suppos_axis_xyz(const O &o, const N &n) {
    MatX<NumType> m = MatX<NumType>::Identity(3, 3);
    double r, x, y, r1, x1, y1, r2, x2, y2, c, s, c1, c2, s1, s2;
    r = std::sqrt(o[0] * o[0] + o[1] * o[1]); x = o[0]; y = o[1];
    if (r != 0) { c = y / r; s = x / r; if (s != 0) m *= z_rot_mat(c, s); }
    r1 = std::sqrt(o[0] * o[0] + o[1] * o[1] + o[2] * o[2]); x1 = std::sqrt(o[0] * o[0] + o[1] * o[1]); y1 = o[2];
    r2 = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]); x2 = std::sqrt(n[0] * n[0] + n[1] * n[1]); y2 = n[2];
    if (r1 != 0 && r2 != 0) {
        c1 = y1 / r1; s1 = x1 / r1; c2 = y2 / r2; s2 = x2 / r2; c = c1*c2 + s1*s2; s = s1*c2 - s2*c1;
        if (s != 0) m *= x_rot_mat(c, s);
    }
    r = std::sqrt(n[0] * n[0] + n[1] * n[1]); x = n[0]; y = n[1];
    if (r != 0) { c = y / r; s = -x / r; if (s != 0) m *= z_rot_mat(c, s); }
    return m;
}

} // namespace geom

} // namespace jnc
