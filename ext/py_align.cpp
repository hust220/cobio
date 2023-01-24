#include "./py_align.h"
#include "./jnc_bio_pdb.h"
#include "./jnc_geom.h"

std::vector<std::array<double, 3>> atoms_coords(PyObject *atoms) {
  int n = PyList_Size(atoms);
  std::vector<std::array<double, 3>> crds(n);
  for (int i = 0; i < n; i++) {
    auto atom = PyList_GET_ITEM(atoms, i);
    auto coord = attr(atom, "coord");
    for (int j = 0; j < 3; j++) {
      crds[i][j] = at(coord, j);
    }
  }
  return crds;
}

std::vector<std::set<std::size_t>> atoms_edges(PyObject *atoms) {
  std::size_t n = PyList_Size(atoms);
  std::vector<std::set<std::size_t>> edges(n);
  auto &&crds = atoms_coords(atoms);
  double cutoff = 1.65 * 1.65;
  for (std::size_t i = 0; i < n; i++) {
    for (std::size_t j = i + 1; j < n; j++) {
      double d2 = jnc::geom::dist2(crds[i], crds[j]);
      if (d2 < cutoff) {
        edges[i].insert(j);
        edges[j].insert(i);
      }
    }
  }
  return edges;
}

struct AtomTypeIdentifier {
  std::vector<std::string> elements{
      "h",  "he", "li", "be", "b",  "c",  "n",  "o",  "f",  "ne", "na", "mg", "al", "si", "p",  "s",  "cl",
      "ar", "k",  "ca", "sc", "ti", "v",  "cr", "mn", "fe", "co", "ni", "cu", "zn", "ga", "ge", "as", "se",
      "br", "kr", "rb", "sr", "y",  "zr", "nb", "mo", "tc", "ru", "rh", "pd", "ag", "cd", "in", "sn", "sb",
      "te", "i",  "xe", "cs", "ba", "la", "ce", "pr", "nd", "pm", "sm", "eu", "gd", "tb", "dy", "ho", "er",
      "tm", "yb", "lu", "hf", "ta", "w",  "re", "os", "ir", "pt", "au", "hg", "tl", "pb", "bi", "po", "at",
      "rn", "fr", "ra", "ac", "th", "pa", "u",  "np", "pu", "am", "cm", "bk", "cf", "es", "fm", "md", "no",
      "lr", "rf", "db", "sg", "bh", "hs", "mt", "ds", "rg", "cn", "nh", "fi", "mc", "lv", "ts", "og"};

  std::map<char, std::map<std::string, int>> ids;

  AtomTypeIdentifier() {
    int n = int(elements.size());
    for (int i = 0; i < n; i++) {
      ids[elements[i][0]][elements[i]] = i;
    }
  }

  int operator()(const std::string &name) {
    int type = -1;
    auto name_l = jnc::string_lower_c(name);
    for (auto &&p : ids[name_l[0]]) {
      auto &key = p.first;
      int nk = int(key.size());
      if (name_l.size() == 1 && nk == 1 && name_l == key) {
        return p.second;
      } else if (name_l.size() > 1) {
        if (nk == 2 && key[0] == name_l[0] && key[1] == name_l[1]) {
          return p.second;
        } else if (nk == 1 && key[0] == name_l[0]) {
          type = p.second;
        }
      }
    }
    return type;
  }
};

jnc::geom::SupPos<double> pdb_align(jnc::bio::Pdb pdb1, jnc::bio::Pdb pdb2) { return jnc::geom::SupPos<double>{}; }

void align_apply(jnc::geom::SupPos<double>, jnc::bio::Pdb pdb) {}

std::vector<int> map_atoms_(PyObject *atoms1, PyObject *atoms2) {
  int n1 = PyList_Size(atoms1);
  int n2 = PyList_Size(atoms2);

  std::vector<std::size_t> nodes1(n1);
  std::vector<std::size_t> nodes2(n2);

  AtomTypeIdentifier ident;
  for (int i = 0; i < n1; i++) {
    auto atom = PyList_GET_ITEM(atoms1, i);
    auto name = o2s(attr(atom, "name"));
    nodes1[i] = std::size_t(ident(name));
  }
  for (int i = 0; i < n2; i++) {
    auto atom = PyList_GET_ITEM(atoms2, i);
    auto name = o2s(attr(atom, "name"));
    nodes2[i] = std::size_t(ident(name));
  }
  auto edges1 = atoms_edges(atoms1);
  auto edges2 = atoms_edges(atoms2);
  std::hash<std::size_t> int_hash;
  auto update_nodes = [](std::vector<std::size_t> &nodes, const std::vector<std::set<std::size_t>> &edges) {
    std::size_t n = nodes.size();
    std::vector<std::size_t> new_nodes(n);
    for (std::size_t i = 0; i < n; i++) {
      std::set<std::size_t> s;
      for (auto &&j : edges[i])
        s.insert(nodes[j]);
      std::size_t hash1 = jnc::hash_range(s.begin(), s.end());
      std::vector<std::size_t> v{nodes[i], hash1};
      new_nodes[i] = jnc::hash_range(v.begin(), v.end());
    }
    for (std::size_t i = 0; i < n; i++)
      nodes[i] = new_nodes[i];
  };
  while (true) {
    //        for (int i = 0; i < n1; i++) { std::cout << nodes1[i] << ' '; }
    //        std::cout << std::endl; for (int i = 0; i < n2; i++) { std::cout
    //        << nodes2[i] << ' '; } std::cout << std::endl; std::cout <<
    //        std::endl;

    std::vector<std::set<std::size_t>> m(n1);
    for (int i = 0; i < n1; i++) {
      update_nodes(nodes1, edges1);
      update_nodes(nodes2, edges2);
    }
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
        if (nodes1[i] == nodes2[j]) {
          m[i].insert(j);
        }
      }
    }

    //        for (int i = 0; i < n1; i++) {
    //            for (auto && j : m[i]) std::cout << j << ' '; std::cout <<
    //            std::endl;
    //        }
    //        std::cout << std::endl;

    auto it = std::find_if(m.begin(), m.end(), [](const std::set<std::size_t> &s) { return s.size() == 0; });
    if (it != m.end()) {
      return std::vector<int>{};
    }

    it = std::find_if(m.begin(), m.end(), [](const std::set<std::size_t> &s) { return s.size() > 1; });
    if (it == m.end()) {
      std::vector<int> v(n1);
      std::transform(m.begin(), m.end(), v.begin(), [](const std::set<std::size_t> &s) { return int(*s.begin()); });
      return v;
    } else {
      std::size_t i = std::distance(m.begin(), it);
      nodes1[i] = int_hash(nodes1[i]);
      std::size_t j = *(m[i].begin());
      nodes2[j] = int_hash(nodes2[j]);
    }
  }
}

int res_natoms(PyObject *res) { return PyList_Size(attr(res, "atoms")); }

PyObject *align(PyObject *self, PyObject *args) {
  const char *seq1, *seq2;
  PARSE_TUPLE(args, "ss", &seq1, &seq2);
  Alignment alignment;
  int score;
  std::string seq1_aligned, seq2_aligned;
  //    std::tie(score, seq1_aligned, seq2_aligned) = alignment(o2s(seq1),
  //    o2s(seq2));
  std::tie(score, seq1_aligned, seq2_aligned) = alignment(seq1, seq2);
  return Py_BuildValue("(N,N,N)", str2obj(seq1_aligned.c_str()), str2obj(seq2_aligned.c_str()), PyLong_FromLong(score));
}

std::map<std::string, std::array<double, 3>> atom_crd_map(PyObject *atoms) {
  int natoms = PyList_Size(atoms);
  std::map<std::string, std::array<double, 3>> rc;
  for (int i = 0; i < natoms; i++) {
    auto atom = PyList_GET_ITEM(atoms, i);
    auto name = o2s(attr(atom, "name"));
    auto coord = attr(atom, "coord");
    std::array<double, 3> crd;
    for (int j = 0; j < 3; j++) {
      crd[j] = at(coord, j);
    }
    rc[name] = crd;
  }
  return std::move(rc);
}

jnc::geom::SupPos<double> _suppos(PyObject *rs1, PyObject *rs2, bool use_common_atoms) {
  int n1 = PyList_Size(rs1);
  int n2 = PyList_Size(rs2);

  jnc::geom::Matd crd1, crd2;

  if (!use_common_atoms) {
    int natoms = 0;
    for (int i = 0; i < n1; i++) {
      auto res = PyList_GET_ITEM(rs1, i);
      natoms += res_natoms(res);
    }

    crd1.resize(natoms, 3);
    crd2.resize(natoms, 3);
    for (int i = 0; i < n1; i++) {
      auto r1 = PyList_GET_ITEM(rs1, i);
      auto r2 = PyList_GET_ITEM(rs2, i);

      auto atoms1 = attr(r1, "atoms");
      auto atoms2 = attr(r2, "atoms");

      int natoms1 = PyList_Size(atoms1);
      for (int i = 0; i < natoms1; i++) {
        auto atom1 = PyList_GET_ITEM(atoms1, i);
        auto atom2 = PyList_GET_ITEM(atoms2, i);
        auto coord1 = attr(atom1, "coord");
        auto coord2 = attr(atom2, "coord");
        for (int j = 0; j < 3; j++) {
          crd1(i, j) = at(coord1, j);
          crd2(i, j) = at(coord2, j);
        }
      }
    }
  } else {
    std::list<std::array<double, 3>> ls1;
    std::list<std::array<double, 3>> ls2;
    for (int i = 0; i < n1; i++) {
      auto r1 = PyList_GET_ITEM(rs1, i);
      auto r2 = PyList_GET_ITEM(rs2, i);

      auto atoms1 = attr(r1, "atoms");
      auto atoms2 = attr(r2, "atoms");

      auto rc1 = atom_crd_map(atoms1);
      auto rc2 = atom_crd_map(atoms2);

      for (auto &&p : rc1) {
        if (rc2.find(p.first) != rc2.end()) {
          ls1.push_back(p.second);
          ls2.push_back(rc2[p.first]);
        }
      }
    }
    int natoms = ls1.size();
    crd1.resize(natoms, 3);
    crd2.resize(natoms, 3);
    auto it1 = ls1.begin();
    auto it2 = ls2.begin();
    for (int i = 0; i < natoms; i++) {
      for (int j = 0; j < 3; j++) {
        crd1(i, j) = (*it1)[j];
        crd2(i, j) = (*it2)[j];
      }
      it1++;
      it2++;
    }
  }
  return jnc::geom::SupPos<double>(crd1, crd2);
}

PyObject *map_atoms(PyObject *self, PyObject *args) {
  PyObject *atoms1, *atoms2;
  PARSE_TUPLE(args, "OO", &atoms1, &atoms2);
  auto &&v = map_atoms_(atoms1, atoms2);
  std::size_t n = v.size();

  PyObject *m = PyList_New(n);
  for (std::size_t i = 0; i < n; i++) {
    PyList_SetItem(m, i, PyLong_FromLong(v[i]));
  }

  return m;
}

PyObject *rmsd(PyObject *self, PyObject *args) {
  PyObject *rs1, *rs2;
  int align;
  int use_common_atoms;
  PARSE_TUPLE(args, "OOpp", &rs1, &rs2, &align, &use_common_atoms);
  double rmsd = 0;
  if (align) {
    auto sp = _suppos(rs1, rs2, use_common_atoms);
    rmsd = sp.rmsd;
  } else {
    int na = 0;
    int n = PyList_Size(rs1);
    for (int i = 0; i < n; i++) {
      auto r1 = PyList_GET_ITEM(rs1, i);
      auto r2 = PyList_GET_ITEM(rs2, i);
      auto atoms1 = attr(r1, "atoms");
      auto atoms2 = attr(r2, "atoms");
      int natoms = PyList_Size(atoms1);
      na += natoms;

      // Map atoms1 to atoms2
      std::vector<int> m(natoms);
      std::iota(m.begin(), m.end(), 0);
      if (use_common_atoms)
        m = map_atoms_(atoms1, atoms2);

      for (int j = 0; j < natoms; j++) {
        auto atom1 = PyList_GET_ITEM(atoms1, j);
        auto atom2 = PyList_GET_ITEM(atoms2, m[j]);
        auto coord1 = attr(atom1, "coord");
        auto coord2 = attr(atom2, "coord");
        for (int k = 0; k < 3; k++) {
          double d = at(coord1, k) - at(coord2, k);
          rmsd += d * d;
        }
      }
    }
    rmsd = std::sqrt(rmsd / na);
  }
  return PyFloat_FromDouble(rmsd);
}

PyObject *suppos(PyObject *self, PyObject *args) {
  PyObject *rs1, *rs2;
  int use_common_atoms;
  PARSE_TUPLE(args, "OOp", &rs1, &rs2, &use_common_atoms);
  auto sp = _suppos(rs1, rs2, use_common_atoms);
  PyObject *rot = PyList_New(3);
  PyObject *c1 = PyList_New(3);
  PyObject *c2 = PyList_New(3);
  for (int i = 0; i < 3; i++) {
    PyList_SetItem(c1, i, PyFloat_FromDouble(sp.c1(i)));
    PyList_SetItem(c2, i, PyFloat_FromDouble(sp.c2(i)));
    PyObject *row = PyList_New(3);
    for (int j = 0; j < 3; j++) {
      PyList_SetItem(row, j, PyFloat_FromDouble(sp.rot(i, j)));
    }
    PyList_SetItem(rot, i, row);
  }
  return Py_BuildValue("((N,N,N),d)", rot, c1, c2, sp.rmsd);
}

// Apply superposition to a list of residues
PyObject *sp_apply(PyObject *self, PyObject *args) {
  PyObject *sp, *rs;
  PARSE_TUPLE(args, "OO", &sp, &rs);

  auto rot = PyTuple_GET_ITEM(sp, 0);
  auto c1 = PyTuple_GET_ITEM(sp, 1);
  auto c2 = PyTuple_GET_ITEM(sp, 2);
  auto rotv = MatViewer(rot);

  int n = PyList_Size(rs);
  for (int i = 0; i < n; i++) {
    auto r = PyList_GET_ITEM(rs, i);
    auto atoms = attr(r, "atoms");
    int natoms = PyList_Size(atoms);
    for (int j = 0; j < natoms; j++) {
      auto atom = PyList_GET_ITEM(atoms, j);
      auto crd = attr(atom, "coord");
      for (int k = 0; k < 3; k++)
        assign(crd, k, at(crd, k) + at(c1, k));
      rotate(crd, rotv);
      for (int k = 0; k < 3; k++)
        assign(crd, k, at(crd, k) + at(c2, k));
    }
  }
  return Py_None;
}

