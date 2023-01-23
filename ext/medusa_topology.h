#pragma once

#include <string>
#include <vector>

namespace medusa {

struct TopAtomType {
  int id;
  double mass;
  std::string name;
};

struct TopAtom {
  int type;
  double charge;
  // std::array<double, 3> ic;
  std::string name;
};

struct TopResidue {
  std::string name;
  std::vector<TopAtom> atoms;
  std::vector<std::array<int, 2>> bonds;
  std::vector<std::array<int, 4>> dihedrals;
  std::vector<std::array<int, 4>> impropers;
};

struct Topology {
  std::string filename;
  std::vector<TopAtomType> atom_table;
  std::vector<TopResidue> resi_table;

  void init(const std::string &filename);
};

} // namespace medusa
