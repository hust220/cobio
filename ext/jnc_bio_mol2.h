#pragma once

#include "./jnc_core.h"

namespace jnc {
namespace bio {

#define MOL2_MOL_TYPE_MAP                                                      \
  ELT(SMALL, /*MOL_TYPE_*/ SMALL)                                              \
  SEP ELT(BIOPOLYMER, /*MOL_TYPE_*/ BIOPOLYMER)                                \
  SEP ELT(PROTEIN, /*MOL_TYPE_*/ PROTEIN)                                      \
  SEP ELT(NUCLEIC_ACID, /*MOL_TYPE_*/ NUCLEIC_ACID)                            \
  SEP ELT(SACCHARIDE, /*MOL_TYPE_*/ SACCHARIDE) SEP

#define MOL2_CHARGE_TYPE_MAP                                                   \
  ELT(NO_CHARGES, /*CHARGE_TYPE_*/ NO_CHARGES)                                 \
  SEP ELT(DEL_RE, /*CHARGE_TYPE_*/ DEL_RE)                                     \
  SEP ELT(GASTEIGER, /*CHARGE_TYPE_*/ GASTEIGER)                               \
  SEP ELT(GAST_HUCK, /*CHARGE_TYPE_*/ GAST_HUCK)                               \
  SEP ELT(HUCKEL, /*CHARGE_TYPE_*/ HUCKEL)                                     \
  SEP ELT(PULLMAN, /*CHARGE_TYPE_*/ PULLMAN)                                   \
  SEP ELT(GAUSS80_CHARGES, /*CHARGE_TYPE_*/ GAUSS80_CHARGES)                   \
  SEP ELT(AMPAC_CHARGES, /*CHARGE_TYPE_*/ AMPAC_CHARGES)                       \
  SEP ELT(MULLIKEN_CHARGES, /*CHARGE_TYPE_*/ MULLIKEN_CHARGES)                 \
  SEP ELT(DICT_CHARGES, /*CHARGE_TYPE_*/ DICT_CHARGES)                         \
  SEP ELT(MMFF94_CHARGES, /*CHARGE_TYPE_*/ MMFF94_CHARGES)                     \
  SEP ELT(USER_CHARGES, /*CHARGE_TYPE_*/ USER_CHARGES) SEP

#define MOL2_BOND_TYPE_MAP                                                     \
  ELT("1", /*BOND_TYPE_*/ SINGLE, single)                                      \
  SEP ELT("2", /*BOND_TYPE_*/ DOUBLE, double)                                  \
      SEP ELT("3", /*BOND_TYPE_*/ TRIPLE, triple)                              \
  SEP ELT("am", /*BOND_TYPE_*/ AMIDE, amide)                                   \
  SEP ELT("ar", /*BOND_TYPE_*/ AROMATIC, aromatic)                             \
  SEP ELT("du", /*BOND_TYPE_*/ DUMMY, dummy)                                   \
  SEP ELT("un", /*BOND_TYPE_*/ UNKNOWN, unknown)                               \
  SEP ELT("nc", /*BOND_TYPE_*/ NOTCONNECTED, not connected) SEP

#define MOL2_BOND_STATUS_MAP                                                   \
  ELT(1, /*BOND_STATUS_*/ TYPECOL) SEP ELT(2, /*BOND_STATUS_*/ GROUP)          \
  SEP ELT(4, /*BOND_STATUS_*/ CAP)                                             \
  SEP ELT(8, /*BOND_STATUS_*/ BACKBONE)                                        \
  SEP ELT(16, /*BOND_STATUS_*/ DICT)                                           \
  SEP ELT(32, /*BOND_STATUS_*/ INTERRES) SEP

#define MOL2_ATOM_TYPE_MAP                                                     \
  ELT(Any, /*ATOM_TYPE_*/ ANY, "any atom")                                     \
  SEP ELT(Hev, /*ATOM_TYPE_*/ HEV, "heavy atom (non hydrogen)") SEP ELT(       \
      Het, /*ATOM_TYPE_*/ HET, "heteroatom = N, O, S, P") SEP                  \
  ELT(LP, /*ATOM_TYPE_*/ LP, "lone pair") SEP                                  \
  ELT(Du, /*ATOM_TYPE_*/ DU, "dummy atom") SEP                                 \
  ELT(Du.C, /*ATOM_TYPE_*/ DU_C, "dummy carbon") SEP                           \
                                                                               \
  ELT(H, /*ATOM_TYPE_*/ H, "hydrogen") SEP                                     \
  ELT(H.spc, /*ATOM_TYPE_*/ H_SPC,                                             \
      "hydrogen in Single Point Charge (SPC) water model") SEP                 \
  ELT(H.t3p, /*ATOM_TYPE_*/ H_T3P,                                             \
      "hydrogen in Transferable intermolecular Potential (TIP3P) water model") \
      SEP                                                                      \
                                                                               \
      ELT(C .3, /*ATOM_TYPE_*/ C_3, "carbon sp3") SEP                          \
      ELT(C .2, /*ATOM_TYPE_*/ C_2, "carbon sp2") SEP                          \
      ELT(C .1, /*ATOM_TYPE_*/ C_1, "carbon sp") SEP                           \
      ELT(C.ar, /*ATOM_TYPE_*/ C_AR, "carbon aromatic") SEP                    \
      ELT(C.cat, /*ATOM_TYPE_*/ C_CAT,                                         \
          "carbocation (C +) used only in a guadinium group") SEP              \
                                                                               \
      ELT(N .4, /*ATOM_TYPE_*/ N_4, "nitrogen sp3 positively charged") SEP     \
      ELT(N .3, /*ATOM_TYPE_*/ N_3, "nitrogen sp3") SEP                        \
      ELT(N .2, /*ATOM_TYPE_*/ N_2, "nitrogen sp2 ") SEP                       \
      ELT(N .1, /*ATOM_TYPE_*/ N_1, "nitrogen sp") SEP                         \
      ELT(N.ar, /*ATOM_TYPE_*/ N_AR, "nitrogen aromatic") SEP                  \
      ELT(N.am, /*ATOM_TYPE_*/ N_AM, "nitrogen amide") SEP                     \
      ELT(N.pl3, /*ATOM_TYPE_*/ N_PL3, "nitrogen trigonal planar") SEP         \
                                                                               \
      ELT(O .3, /*ATOM_TYPE_*/ O_3, "oxygen sp3") SEP                          \
      ELT(O .2, /*ATOM_TYPE_*/ O_2, "oxygen sp2") SEP                          \
      ELT(O.co2, /*ATOM_TYPE_*/ O_CO2,                                         \
          "oxygen in carboxylate and phosphate groups") SEP                    \
      ELT(O.spc, /*ATOM_TYPE_*/ O_SPC,                                         \
          "oxygen in Single Point Charge (SPC) water model ") SEP              \
      ELT(O.t3p, /*ATOM_TYPE_*/ O_T3P,                                         \
          "oxygen in Transferable Intermolecular Potential (TIP3P) water "     \
          "model") SEP                                                         \
                                                                               \
      ELT(S .3, /*ATOM_TYPE_*/ S_3, "sulfur sp3") SEP                          \
      ELT(S .2, /*ATOM_TYPE_*/ S_2, "sulfur sp2") SEP                          \
      ELT(S.O, /*ATOM_TYPE_*/ S_O, "sulfoxide sulfur") SEP                     \
      ELT(S.O2, /*ATOM_TYPE_*/ S_O2, "sulfone sulfur") SEP                     \
      ELT(P .3, /*ATOM_TYPE_*/ P_3, "phosphorous sp3") SEP                     \
                                                                               \
      ELT(Hal, /*ATOM_TYPE_*/ HAL,                                             \
          "halogen: fluorine(F), chlorine(Cl), bromine(Br), iodine(I), "       \
          "astatine(At)") SEP                                                  \
      ELT(F, /*ATOM_TYPE_*/ F, "fluorine") SEP                                 \
      ELT(Cl, /*ATOM_TYPE_*/ CL, "chlorine") SEP                               \
      ELT(Br, /*ATOM_TYPE_*/ BR, "bromine") SEP                                \
      ELT(I, /*ATOM_TYPE_*/ I, "iodine") SEP                                   \
                                                                               \
      ELT(Li, /*ATOM_TYPE_*/ LI, "lithium") SEP                                \
      ELT(Na, /*ATOM_TYPE_*/ NA, "sodium") SEP                                 \
      ELT(Mg, /*ATOM_TYPE_*/ MG, "magnesium") SEP                              \
      ELT(Fe, /*ATOM_TYPE_*/ FE, "iron") SEP                                   \
      ELT(Al, /*ATOM_TYPE_*/ AL, "aluminum") SEP                               \
      ELT(Si, /*ATOM_TYPE_*/ SI, "silicon") SEP                                \
      ELT(K, /*ATOM_TYPE_*/ K, "potassium") SEP                                \
      ELT(Ca, /*ATOM_TYPE_*/ CA, "calcium") SEP                                \
      ELT(Cr.th, /*ATOM_TYPE_*/ CR_TH, "chromium (tetrahedral)") SEP           \
      ELT(Cr.oh, /*ATOM_TYPE_*/ CR_OH, "chromium (octahedral)") SEP            \
      ELT(Co.oh, /*ATOM_TYPE_*/ CO_OH, "cobalt (octahedral)") SEP              \
      ELT(Mn, /*ATOM_TYPE_*/ MN, "manganese") SEP                              \
      ELT(Cu, /*ATOM_TYPE_*/ CU, "copper") SEP                                 \
      ELT(Zn, /*ATOM_TYPE_*/ ZN, "zinc") SEP                                   \
      ELT(Se, /*ATOM_TYPE_*/ SE, "selenium") SEP                               \
      ELT(Mo, /*ATOM_TYPE_*/ MO, "molybdenum") SEP                             \
      ELT(Sn, /*ATOM_TYPE_*/ SN, "tin") SEP

#define MOL2_ATOM_STATUS_MAP                                                   \
  ELT(1, /*ATOM_STATUS_*/ DSPMOD) SEP ELT(2, /*ATOM_STATUS_*/ TYPECOL)         \
  SEP ELT(4, /*ATOM_STATUS_*/ CAP)                                             \
  SEP ELT(8, /*ATOM_STATUS_*/ BACKBONE)                                        \
  SEP ELT(16, /*ATOM_STATUS_*/ DICT)                                           \
  SEP ELT(32, /*ATOM_STATUS_*/ ESSENTIAL)                                      \
  SEP ELT(64, /*ATOM_STATUS_*/ WATER)                                          \
  SEP ELT(128, /*ATOM_STATUS_*/ DIRECT) SEP

#define JN_MOL2_ATOM_TABLE                                                     \
  ELT("H", ATOM_TYPE_H, ATOM_ELEMENT_H)                                        \
  SEP ELT("H.SPC", ATOM_TYPE_H_SPC, ATOM_ELEMENT_H)                            \
  SEP ELT("H.T3P", ATOM_TYPE_H_T3P, ATOM_ELEMENT_H)                            \
  SEP                                                                          \
                                                                               \
      ELT("C.1", ATOM_TYPE_C_1, ATOM_ELEMENT_C)                                \
  SEP ELT("C.2", ATOM_TYPE_C_2, ATOM_ELEMENT_C)                                \
  SEP ELT("C.3", ATOM_TYPE_C_3, ATOM_ELEMENT_C)                                \
  SEP ELT("C.AR", ATOM_TYPE_C_AR, ATOM_ELEMENT_C)                              \
  SEP ELT("C.CAT", ATOM_TYPE_C_CAT, ATOM_ELEMENT_C)                            \
  SEP                                                                          \
                                                                               \
      ELT("N.1", ATOM_TYPE_N_1, ATOM_ELEMENT_N)                                \
  SEP ELT("N.2", ATOM_TYPE_N_2, ATOM_ELEMENT_N)                                \
  SEP ELT("N.3", ATOM_TYPE_N_3, ATOM_ELEMENT_N)                                \
  SEP ELT("N.4", ATOM_TYPE_N_4, ATOM_ELEMENT_N)                                \
  SEP ELT("N.AR", ATOM_TYPE_N_AR, ATOM_ELEMENT_N)                              \
  SEP ELT("N.AM", ATOM_TYPE_N_AM, ATOM_ELEMENT_N)                              \
  SEP ELT("N.PL3", ATOM_TYPE_N_PL3, ATOM_ELEMENT_N)                            \
  SEP                                                                          \
                                                                               \
      ELT("O.2", ATOM_TYPE_O_2, ATOM_ELEMENT_O)                                \
  SEP ELT("O.3", ATOM_TYPE_O_3, ATOM_ELEMENT_O)                                \
  SEP ELT("O.CO2", ATOM_TYPE_O_CO2, ATOM_ELEMENT_O)                            \
  SEP ELT("O.T3P", ATOM_TYPE_O_T3P, ATOM_ELEMENT_O)                            \
  SEP ELT("O.SPC", ATOM_TYPE_O_SPC, ATOM_ELEMENT_O)                            \
  SEP                                                                          \
                                                                               \
      ELT("P.3", ATOM_TYPE_P_3, ATOM_ELEMENT_P)                                \
  SEP                                                                          \
                                                                               \
      ELT("S.2", ATOM_TYPE_S_2, ATOM_ELEMENT_S)                                \
  SEP ELT("S.3", ATOM_TYPE_S_3, ATOM_ELEMENT_S)                                \
  SEP ELT("S.O", ATOM_TYPE_S_O, ATOM_ELEMENT_S)                                \
  SEP ELT("S.O2", ATOM_TYPE_S_O2, ATOM_ELEMENT_S)                              \
  SEP                                                                          \
                                                                               \
      ELT("K", ATOM.TYPE_K, ATOM_ELEMENT_K)                                    \
  SEP ELT("F", ATOM.TYPE_F, ATOM_ELEMENT_F)                                    \
  SEP ELT("I", ATOM.TYPE_I, ATOM_ELEMENT_I)                                    \
  SEP ELT("LP", ATOM.TYPE_LP, ATOM_ELEMENT_LP)                                 \
  SEP ELT("MG", ATOM.TYPE_MG, ATOM_ELEMENT_MG)                                 \
  SEP ELT("CR.OH", ATOM_TYPE_CR_OH, ATOM_ELEMENT_CR)                           \
  SEP ELT("CR.TH", ATOM_TYPE_CR_TH, ATOM_ELEMENT_CR)                           \
  SEP ELT("CO.OH", ATOM_TYPE_CO_OH, ATOM_ELEMENT_CO)                           \
  SEP ELT("CL", ATOM.TYPE_CL, ATOM_ELEMENT_CL)                                 \
  SEP ELT("SE", ATOM.TYPE_SE, ATOM_ELEMENT_SE)                                 \
  SEP ELT("NA", ATOM.TYPE_NA, ATOM_ELEMENT_NA)                                 \
  SEP ELT("SI", ATOM.TYPE_SI, ATOM_ELEMENT_SI)                                 \
  SEP ELT("FE", ATOM.TYPE_FE, ATOM_ELEMENT_FE)                                 \
  SEP ELT("ZN", ATOM.TYPE_ZN, ATOM_ELEMENT_ZN)                                 \
  SEP ELT("SN", ATOM.TYPE_SN, ATOM_ELEMENT_SN)                                 \
  SEP ELT("LI", ATOM.TYPE_LI, ATOM_ELEMENT_LI)                                 \
  SEP ELT("AL", ATOM.TYPE_AL, ATOM_ELEMENT_AL)                                 \
  SEP ELT("CA", ATOM.TYPE_CA, ATOM_ELEMENT_CA)                                 \
  SEP ELT("MN", ATOM.TYPE_MN, ATOM_ELEMENT_MN)                                 \
  SEP ELT("CU", ATOM.TYPE_CU, ATOM_ELEMENT_CU)                                 \
  SEP ELT("BR", ATOM.TYPE_BR, ATOM_ELEMENT_BR)                                 \
  SEP ELT("MO", ATOM.TYPE_MO, ATOM_ELEMENT_MO)                                 \
  SEP                                                                          \
                                                                               \
      ELT("DU", ATOM.TYPE_DU, ATOM_ELEMENT_X)                                  \
  SEP ELT("ANY", ATOM.TYPE_ANY, ATOM_ELEMENT_X)                                \
  SEP ELT("DU.C", ATOM_TYPE_DU_C, ATOM_ELEMENT_X)                              \
  SEP ELT("HAL", ATOM.TYPE_HAL, ATOM_ELEMENT_X)                                \
  SEP ELT("HET", ATOM.TYPE_HET, ATOM_ELEMENT_X)                                \
  SEP ELT("HEV", ATOM.TYPE_HEV, ATOM_ELEMENT_X)

/**
 * Mol2MolType
 */
class Mol2MolType {
public:
  enum {
#define ELT(name, type) type
#define SEP ,
    MOL2_MOL_TYPE_MAP
#undef ELT
#undef SEP
  };

  static std::string get_name(int id) {
    std::map<int, std::string> m{
#define ELT(name, type)                                                        \
  { type, #name }
#define SEP ,
        MOL2_MOL_TYPE_MAP
#undef ELT
#undef SEP
    };
    return m[id];
  }

  static int get_id(const std::string &name) {
    static std::map<std::string, int> m{
#define ELT(name, type)                                                        \
  { #name, type }
#define SEP ,
        MOL2_MOL_TYPE_MAP
#undef ELT
#undef SEP
    };
    return m[name];
  }
};

/**
 * Mol2 Charge Type
 */
class Mol2ChargeType {
public:
  enum {
#define ELT(name, type) type
#define SEP ,
    MOL2_CHARGE_TYPE_MAP
#undef ELT
#undef SEP
  };

  static std::string get_name(int id) {
    std::map<int, std::string> m{
#define ELT(name, type)                                                        \
  { type, #name }
#define SEP ,
        MOL2_CHARGE_TYPE_MAP
#undef ELT
#undef SEP
    };
    return m[id];
  }

  static int get_id(const std::string &name) {
    std::map<std::string, int> m{
#define ELT(name, type)                                                        \
  { #name, type }
#define SEP ,
        MOL2_CHARGE_TYPE_MAP
#undef ELT
#undef SEP
    };
    return m[name];
  }
};

/**
 * Mol2 Bond Type
 */
class Mol2BondType {
public:
  enum {
#define ELT(name, type, annotation) type
#define SEP ,
    MOL2_BOND_TYPE_MAP
#undef ELT
#undef SEP
  };

  static std::string get_name(int id) {
    std::map<int, std::string> m{
#define ELT(name, type, annotation)                                            \
  { type, name }
#define SEP ,
        MOL2_BOND_TYPE_MAP
#undef ELT
#undef SEP
    };
    return m[id];
  }

  static int get_id(const std::string &name) {
    std::map<std::string, int> m{
#define ELT(name, type, annotation)                                            \
  { name, type }
#define SEP ,
        MOL2_BOND_TYPE_MAP
#undef ELT
#undef SEP
    };
    return m[name];
  }
};

/**
 * Mol2BondStatus
 */
class Mol2BondStatus {
public:
  enum {
#define ELT(n, status) status = n
#define SEP ,
    MOL2_BOND_STATUS_MAP
#undef ELT
#undef SEP
  };

  static std::string get_name(int id) {
    std::map<int, std::string> m{
#define ELT(n, status)                                                         \
  { status, #status }
#define SEP ,
        MOL2_BOND_STATUS_MAP
#undef ELT
#undef SEP
    };
    return m[id];
  }

  static int get_id(const std::string &name) {
    std::map<std::string, int> m{
#define ELT(n, status)                                                         \
  { #status, status }
#define SEP ,
        MOL2_BOND_STATUS_MAP
#undef ELT
#undef SEP
    };
    return m[name];
  }
};

/**
 * Mol2AtomType
 */
class Mol2AtomType {
public:
  enum {
#define ELT(name, type, annotation) type
#define SEP ,
    MOL2_ATOM_TYPE_MAP
#undef ELT
#undef SEP
  };

  static std::string get_name(int id) {
    std::map<int, std::string> m{
#define ELT(name, type, annotation)                                            \
  { type, #name }
#define SEP ,
        MOL2_ATOM_TYPE_MAP
#undef ELT
#undef SEP
    };
    return m[id];
  }

  static int get_id(const std::string &name) {
    std::map<std::string, int> m{
#define ELT(name, type, annotation)                                            \
  { #name, type }
#define SEP ,
        MOL2_ATOM_TYPE_MAP
#undef ELT
#undef SEP
    };
    return m[name];
  }
};

/**
 * Mol2AtomStatus related definitions
 */
class Mol2AtomStatus {
public:
  enum {
#define ELT(n, status) status = n
#define SEP ,
    MOL2_ATOM_STATUS_MAP
#undef ELT
#undef SEP
  };

  static std::string get_name(int id) {
    std::map<int, std::string> m{
#define ELT(n, status)                                                         \
  { status, #status }
#define SEP ,
        MOL2_ATOM_STATUS_MAP
#undef ELT
#undef SEP
    };
    return m[id];
  }

  static int get_id(const std::string &name) {
    std::map<std::string, int> m{
#define ELT(n, status)                                                         \
  { #status, status }
#define SEP ,
        MOL2_ATOM_STATUS_MAP
#undef ELT
#undef SEP
    };
    return m[name];
  }
};

struct Mol2Bond {
  int origin_atom_id;
  int target_atom_id;
  std::string type;
  std::string status_bits = "";
};
using Mol2Bonds = Vec<Mol2Bond>;

struct Mol2Atom {
  std::string name;
  double x;
  double y;
  double z;
  std::string type;
  int subst_id = INT_MAX;
  std::string subst_name = "";
  double charge = DBL_MAX;
  std::string status_bit = "";

  double &operator[](int i) { return i == 0 ? x : (i == 1 ? y : z); }

  const double &operator[](int i) const {
    return i == 0 ? x : (i == 1 ? y : z);
  }
};
using Mol2Atoms = Vec<Mol2Atom>;

struct Mol2Substructure {
  std::string name;
  int root_atom;
  std::string type;
  int dict_type = INT_MAX;
  std::string chain;
  std::string sub_type;
  int inter_bonds = INT_MAX;
  std::string status;
  std::string comment;
};
using Mol2Substructures = Vec<Mol2Substructure>;

struct Mol2 {
  std::string mol_name;

  int num_atoms;
  int num_bonds = INT_MAX;
  int num_subst = INT_MAX;
  int num_feat = INT_MAX;
  int num_sets = INT_MAX;

  std::string mol_type;
  std::string charge_type;
  std::string status_bits = "";
  std::string mol_comment = "";

  Mol2Atoms atoms;
  Mol2Bonds bonds;
  Mol2Substructures substructures;

  void remove_hydrogens() {
    Mol2Atoms as;
    Vec<int> vi(atoms.size(), -1);
    int i = 0;
    int j = 0;
    for (auto &&atom : atoms) {
      auto &type = atom.type;
      if (type != "H" && type.find("H.") == std::string::npos) {
        vi[i] = j;
        as.push_back(atom);
        j++;
      }
      i++;
    }
    atoms = std::move(as);

    Mol2Bonds bs;
    for (auto &&bond : bonds) {
      bond.origin_atom_id = vi[bond.origin_atom_id];
      bond.target_atom_id = vi[bond.target_atom_id];
      if (bond.origin_atom_id != -1 && bond.target_atom_id != -1)
        bs.push_back(bond);
    }
    bonds = std::move(bs);

    num_atoms = atoms.size();
    num_bonds = bonds.size();
  }

  void write(std::string fn) {
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

    if (!status_bits.empty())
      out << status_bits << std::endl;
    if (!mol_comment.empty())
      out << mol_comment << std::endl;
  }

  void write_atoms(std::ostream &out) {
    out << "@<TRIPOS>ATOM" << std::endl;
    int i = 0;
    for (auto &&atom : atoms) {
      out << i + 1 << ' ' << atom.name << ' ' << atom.x << ' ' << atom.y << ' '
          << atom.z << ' ' << atom.type;
      if (atom.subst_id != INT_MAX) {
        out << ' ' << atom.subst_id + 1;
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
    if (bonds.empty())
      return;

    out << "@<TRIPOS>BOND" << std::endl;
    int i = 0;
    for (auto &&bond : bonds) {
      out << i + 1 << ' ' << bond.origin_atom_id + 1 << ' '
          << bond.target_atom_id + 1 << ' ' << bond.type;
      if (!bond.status_bits.empty())
        out << ' ' << bond.status_bits;
      out << std::endl;
      i++;
    }
  }

  void write_substructures(std::ostream &out) {
    if (substructures.empty())
      return;

    out << "@<TRIPOS>SUBSTRUCTURE" << std::endl;
    int i = 0;
    for (auto &&subst : substructures) {
      out << i + 1 << ' ' << subst.name << ' ' << subst.root_atom;
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

inline Vec<Mol2> read_mol2s(const std::string &fn) {
  std::string section;
  Vec<Mol2> ms;
  Mol2 *p = NULL;

  int mol_lines = 0;
  file_each_line(fn, [&section, &mol_lines, &p, &ms](std::string &line) {
    if (!line.compare(0, 9, "@<TRIPOS>")) {
      section = line.substr(9);
      if (section == "MOLECULE") {
        ms.push_back(Mol2{});
        p = &(ms.back());
      }
    } else if (string_tokenize(line, " \t").empty()) {
    } else {
      if (section == "MOLECULE") {
        mol_lines++;
        if (mol_lines == 1)
          p->mol_name = string_trim_c(line);
        else if (mol_lines == 2) {
          auto &&vs = string_tokenize(line, " \t");
          if (vs.size() >= 1)
            p->num_atoms = JN_INT(vs[0]);
          if (vs.size() >= 2)
            p->num_bonds = JN_INT(vs[1]);
          if (vs.size() >= 3)
            p->num_subst = JN_INT(vs[2]);
          if (vs.size() >= 4)
            p->num_feat = JN_INT(vs[3]);
          if (vs.size() >= 5)
            p->num_sets = JN_INT(vs[4]);
        } else if (mol_lines == 3) {
          p->mol_type = string_trim_c(line);
        } else if (mol_lines == 4)
          p->charge_type = string_trim_c(line);
        else if (mol_lines == 5)
          p->status_bits = string_trim_c(line);
        else if (mol_lines == 6)
          p->mol_comment = string_trim_c(line);
      } else if (section == "ATOM") {
        auto &&vs = string_tokenize(line, " \t");
        Mol2Atom atom;
        if (vs.size() >= 6) {
          atom.name = vs[1];
          atom.x = JN_DBL(vs[2]);
          atom.y = JN_DBL(vs[3]);
          atom.z = JN_DBL(vs[4]);
          atom.type = vs[5];
          if (vs.size() >= 7)
            atom.subst_id = JN_INT(vs[6]) - 1;
          if (vs.size() >= 8)
            atom.subst_name = vs[7];
          if (vs.size() >= 9)
            atom.charge = JN_DBL(vs[8]);
          if (vs.size() >= 10)
            atom.status_bit = vs[9];

          p->atoms.push_back(atom);
        } else {
          JN_DIE("Wrong mol2 format!");
        }
      } else if (section == "BOND") {
        auto &&vs = string_tokenize(line, " \t");
        if (vs.size() >= 4) {
          Mol2Bond bond;
          bond.origin_atom_id = JN_INT(vs[1]) - 1;
          bond.target_atom_id = JN_INT(vs[2]) - 1;
          bond.type = vs[3];
          if (vs.size() > 4) {
            bond.status_bits = vs[4];
          }
          p->bonds.push_back(bond);
        }
      } else if (section == "SUBSTRUCTURE") {
        auto &&vs = string_tokenize(line, " \t");
        if (vs.size() >= 3) {
          Mol2Substructure subst;
          subst.name = vs[1];
          subst.root_atom = JN_INT(vs[2]);
          if (vs.size() > 3)
            subst.type = vs[3];
          if (vs.size() > 4)
            subst.dict_type = JN_INT(vs[4]);
          if (vs.size() > 5)
            subst.chain = vs[5];
          if (vs.size() > 6)
            subst.sub_type = vs[6];
          if (vs.size() > 7)
            subst.inter_bonds = JN_INT(vs[7]);
          if (vs.size() > 8)
            subst.status = vs[8];
          if (vs.size() > 9)
            subst.comment = vs[9];

          p->substructures.push_back(subst);
        }
      }
    }
    return JN_GO;
  });
  return std::move(ms);
}

} // namespace bio
} // namespace jnc
