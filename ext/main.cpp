#define PY_SSIZE_T_CLEAN

#include "Python.h"
#include "jnc.h"
#include "align.h"
#include <numeric>
#include <stdio.h>
#include <array>
#include <vector>
#include <unordered_set>

//aa1 = ['A',   'R',   'N',   'D',   'C',   'Q',   'E',   'G',   'H',   'I',   'L',   'K',   'M',   'F',   'P',   'S',   'T',   'W',   'Y',   'V',   'O',   'U',   'B',   'Z',   'X',   'J',   '*']
//aa2 = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'PYL', 'SEC', 'ASX', 'GLX', 'XAA', 'XLE', 'TERM']

struct MatViewer {
    PyObject *obj;
    MatViewer(PyObject *o) : obj(o) {}
    double operator()(int i, int j) const {
        return PyFloat_AS_DOUBLE(PyList_GET_ITEM(PyList_GET_ITEM(obj, i), j));
    }
};

double at(PyObject *vec, int i) {
    return PyFloat_AS_DOUBLE(PyList_GET_ITEM(vec, i));
}

void assign(PyObject *vec, int i, double v) {
    PyList_SET_ITEM(vec, i, PyFloat_FromDouble(v));
}

char *o2s(PyObject *o) {
    return PyBytes_AS_STRING(PyUnicode_AsUTF8String(o));
}

int o2i(PyObject *o) {
    return int(PyLong_AsLong(PyNumber_Long(o)));
}

double o2d(PyObject *o) {
    return PyFloat_AS_DOUBLE(PyNumber_Float(o));
}

template<typename Array_>
PyObject *a2o(const Array_ &c, int size) {
    PyObject *obj = PyList_New(size);
    for (int i = 0; i < size; i++) {
        PyList_SET_ITEM(obj, i, PyFloat_FromDouble(c[i]));
    }
    return obj;
}

using Residue = jnc::pdb::Residue;
using Atom = jnc::pdb::Atom;

Residue obj2residue(PyObject *o) {
    Residue residue;
    return std::move(residue);
}

PyObject *attr(PyObject *o, const char *key) {
    return PyObject_GetAttrString(o, key);
}

Atom obj2atom(PyObject *o) {
    Atom atom;
    atom.name = o2s(attr(o, "name"));
    auto coord = attr(o, "coord");
    atom[0] = PyFloat_AS_DOUBLE(PyList_GET_ITEM(coord, 0));
    atom[1] = PyFloat_AS_DOUBLE(PyList_GET_ITEM(coord, 1));
    atom[2] = PyFloat_AS_DOUBLE(PyList_GET_ITEM(coord, 2));
    return std::move(atom);
}

#define Err(...) PyErr_SetString(PyExc_TypeError, __VA_ARGS__)

#define PARSE_TUPLE(...) \
    if (!PyArg_ParseTuple(__VA_ARGS__)) { \
        Err("parameter type error"); \
        return NULL; \
    } \

static PyObject *add(PyObject *i, PyObject *j) {
    return PyLong_FromLong(0);
}

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
    return std::move(crds);
}

std::vector<std::set<std::size_t>> atoms_edges(PyObject *atoms) {
    std::size_t n = PyList_Size(atoms);
    std::vector<std::set<std::size_t>> edges(n);
    auto &&crds = atoms_coords(atoms);
    double cutoff = 1.65*1.65;
    for (std::size_t i = 0; i < n; i++) {
        for (std::size_t j = i + 1; j < n; j++) {
            double d2 = jnc::geom::dist2(crds[i], crds[j]);
            if (d2 < cutoff) {
                edges[i].insert(j);
                edges[j].insert(i);
            }
        }
    }
    return std::move(edges);
}

struct AtomTypeIdentifier {
    std::vector<std::string> elements {
        "h", "he", "li", "be", "b", "c", "n", "o", "f", "ne", "na", "mg", "al", "si", "p", "s", "cl", "ar", "k", "ca", "sc", "ti", "v", "cr", "mn", "fe", "co", "ni", "cu", "zn", "ga", "ge", "as", "se", "br", "kr", "rb", "sr", "y", "zr", "nb", "mo", "tc", "ru", "rh", "pd", "ag", "cd", "in", "sn", "sb", "te", "i", "xe", "cs", "ba", "la", "ce", "pr", "nd", "pm", "sm", "eu", "gd", "tb", "dy", "ho", "er", "tm", "yb", "lu", "hf", "ta", "w", "re", "os", "ir", "pt", "au", "hg", "tl", "pb", "bi", "po", "at", "rn", "fr", "ra", "ac", "th", "pa", "u", "np", "pu", "am", "cm", "bk", "cf", "es", "fm", "md", "no", "lr", "rf", "db", "sg", "bh", "hs", "mt", "ds", "rg", "cn", "nh", "fi", "mc", "lv", "ts", "og"
    };
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
        for (auto && p : ids[name_l[0]]) {
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
            for (auto && j : edges[i]) s.insert(nodes[j]);
            std::size_t hash1 = jnc::hash_range(s.begin(), s.end());
            std::vector<std::size_t> v{nodes[i], hash1};
            new_nodes[i] = jnc::hash_range(v.begin(), v.end());
        }
        for (std::size_t i = 0; i < n; i++) nodes[i] = new_nodes[i];
    };
    while (true) {
//        for (int i = 0; i < n1; i++) { std::cout << nodes1[i] << ' '; } std::cout << std::endl;
//        for (int i = 0; i < n2; i++) { std::cout << nodes2[i] << ' '; } std::cout << std::endl;
//        std::cout << std::endl;

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
//            for (auto && j : m[i]) std::cout << j << ' '; std::cout << std::endl;
//        }
//        std::cout << std::endl;

        auto it = std::find_if(m.begin(), m.end(), [](const std::set<std::size_t> &s){
            return s.size() == 0;
        });
        if (it != m.end()) {
            return std::vector<int>{};
        }

        it = std::find_if(m.begin(), m.end(), [](const std::set<std::size_t> &s){
            return s.size() > 1;
        });
        if (it == m.end()) {
            std::vector<int> v(n1);
            std::transform(m.begin(), m.end(), v.begin(), [](const std::set<std::size_t> &s){return int(*s.begin());});
            return v;
        } else {
            std::size_t i = std::distance(m.begin(), it);
            nodes1[i] = int_hash(nodes1[i]);
            std::size_t j = *(m[i].begin());
            nodes2[j] = int_hash(nodes2[j]);
        }
    }
}

std::string model_seq(const jnc::pdb::Model &m) {
    auto rs = m.presidues();
    int n = rs.size();
    std::string seq(n, 'X');
    for (int i = 0; i < n; i++) {
        seq[i] = jnc::string_trim_c(rs[i]->name)[0];
    }
    return seq;
}

jnc::geom::SupPos<double> pdb_align(jnc::pdb::Pdb pdb1, jnc::pdb::Pdb pdb2) {
    return jnc::geom::SupPos<double>{};
}

void align_apply(jnc::geom::SupPos<double>, jnc::pdb::Pdb pdb) {
}

//static PyObject *copy_atom(PyObject *atom) {
//    auto type = atom->ob_type;
//    auto n = type->tp_alloc(type, 0);
//
//    auto && name = attr(atom, "name");
//    auto && c = attr(atom, "coord");
//
//    PyObject_SetAttrString(n, "name", PyUnicode_FromString(name));
//    PyObject *coord = PyList_New(3);
//    for (int i = 0; i < 3; i++) {
//        PyList_SET_ITEM(coord, i, PyLong_FromLong(at(c, i)));
//    }
//    PyObject_SetAttr(n, "coord", coord);
//
//    return n;
//}

//static PyObject *copy_atoms(PyObject *atoms) {
//    int sz = PyList_Size(atoms);
//    PyObject *as = PyList_New(sz);
//    for (int i = 0; i < sz; i++) {
//        PyList_SET_ITEM(as, i, copy_atom(PyList_GET_ITEM(atoms, i)));
//    }
//    return as;
//}

//static PyObject *newatom(PyObject *self, PyObject *args) {
//    PyObject *atom;
//    PARSE_TUPLE(args, "O", &atom);
//    return copy_atom(atom);
//}

int res_natoms(PyObject *res) {
    return PyList_Size(attr(res, "atoms"));
}

PyObject *str2obj(const char *str) {
//    return PyUnicode_DecodeUTF8(PyBytes_FromString(str));
    return Py_BuildValue("s", str);
}

static PyObject *align(PyObject *self, PyObject *args) {
    const char *seq1, *seq2;
    PARSE_TUPLE(args, "ss", &seq1, &seq2);
    Alignment alignment;
    int score;
    std::string seq1_aligned, seq2_aligned;
//    std::tie(score, seq1_aligned, seq2_aligned) = alignment(o2s(seq1), o2s(seq2));
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

static jnc::geom::SupPos<double> _suppos(PyObject *rs1, PyObject *rs2, bool use_common_atoms) {
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

            for (auto && p : rc1) {
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

static PyObject *map_atoms(PyObject *self, PyObject *args) {
    PyObject *atoms1, *atoms2;
    PARSE_TUPLE(args, "OO", &atoms1, &atoms2);
    auto && v = map_atoms_(atoms1, atoms2);
    std::size_t n = v.size();

    PyObject *m = PyList_New(n);
    for (std::size_t i = 0; i < n; i++) {
        PyList_SetItem(m, i, PyLong_FromLong(v[i]));
    }

    return m;
}

static PyObject *rmsd(PyObject *self, PyObject *args) {
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
            if (use_common_atoms) m = map_atoms_(atoms1, atoms2);

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

static PyObject *suppos(PyObject *self, PyObject *args) {
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
            PyList_SetItem(row, j, PyFloat_FromDouble(sp.rot(i,j)));
        }
        PyList_SetItem(rot, i, row);
    }
    return Py_BuildValue("((N,N,N),d)", rot, c1, c2, sp.rmsd);
}


template<typename _Mat>
void rotate(PyObject *o, const _Mat &mat) {
    double _x = at(o, 0);
    double _y = at(o, 1);
    double _z = at(o, 2);
    assign(o, 0, _x * mat(0, 0) + _y * mat(1, 0) + _z * mat(2, 0));
    assign(o, 1, _x * mat(0, 1) + _y * mat(1, 1) + _z * mat(2, 1));
    assign(o, 2, _x * mat(0, 2) + _y * mat(1, 2) + _z * mat(2, 2));
}

// Apply superposition to a list of residues
static PyObject *sp_apply(PyObject *self, PyObject *args) {
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
            for (int k = 0; k < 3; k++) assign(crd, k, at(crd, k)+at(c1, k));
            rotate(crd, rotv);
            for (int k = 0; k < 3; k++) assign(crd, k, at(crd, k)+at(c2, k));
        }
    }
    return Py_None;
}

//////////////////////////////////////////////////////////////////
// Pocket

struct PocketReceptorAtom {
    std::string name;
    std::array<double, 3> coord;
    std::array<int, 3> index; // index in the grid
    int isLigand;
};

struct PocketLigandAtom {
    std::array<int, 3> index;
    uint element = 0; // 0: None, 1: C, 2: 
    PocketLigandAtom *prev = NULL;
};

// static std::string pocket_ligand_atom_element(uint element) {
//     const static std::vector<std::string> m = {"X", "C", "N", "O", "P", "S"};
//     if (element < 0 || element >= m.size()) {
//         Err("Atom element not supported!");
//         return "";
//     } else {
//         return m[element];
//     }
// };

struct PocketGrid {
    std::array<double, 3> center;
    double box;
    double bin;
    int size;
};

template<typename Grid_, typename Coord_>
std::array<int, 3> grid_c2i(const Grid_ &grid, const Coord_ &c) {
    std::array<int, 3> ind;
    for (int i = 0; i < 3; i++) ind[i] = int((c[i] - (grid.center[i] - grid.box / 2.0)) / grid.bin);
    return ind;
}

template<typename Grid_, typename Ind_>
std::array<double, 3> grid_i2c(const Grid_ &grid, const Ind_ &ind) {
    std::array<double, 3> c;
    for (int i = 0; i < 3; i++) c[i] = ind[i] * grid.bin + (grid.center[i] - grid.box / 2.0);
    return c;
}

struct PocketReceptor {
    std::vector<PocketReceptorAtom> atoms;
};

struct PocketLigand {
    PocketLigandAtom lastAtom; // Only the last atom is owned by this ligand.
};

struct Pocket {
    std::size_t id = 0;
    PocketGrid *grid = NULL; // shared by all the nodes
    PocketReceptor *receptor = NULL; // shared by all the nodes
    PocketLigand *ligand = NULL; // owned by the current node
};

template<typename Atom1_, typename Atom2_>
bool atom_in_box(Atom1_ &&atom1, Atom2_ &&atom2, double box) {
    double halfBox = box / 2.0;
    double x = atom1[0] - atom2[0];
    double y = atom1[1] - atom2[1];
    double z = atom1[2] - atom2[2];
    return std::abs(x) < halfBox && std::abs(y) < halfBox && std::abs(z) < halfBox;
}

template<typename Residue_>
std::array<double, 3> residue_center(const Residue_ &residue) {
    std::array<double, 3> c {0, 0, 0};
    int iatom = 0;
    for (auto && atom : residue) {
        for (int i = 0; i < 3; i++) {
            c[i] += atom[i];
        }
        iatom += 1;
    }
    for (int i = 0; i < 3; i++) {
        c[i] /= iatom;
    }
    return c;
}

typedef struct {
    PyObject_HEAD
    Pocket *pocket;
} PocketObject;

static PyObject *Pocket_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    PocketObject *self;
    self = (PocketObject *) type->tp_alloc(type, 0);
    if (self != NULL) {
        self->pocket = new Pocket;
        self->pocket->receptor = NULL;
        self->pocket->receptor = NULL;
        self->pocket->ligand = NULL;
        self->pocket->id = 0;
    }
    return (PyObject *) self;
}

static int Pocket_init(PocketObject *self, PyObject *args, PyObject *kwargs) {
    const char *receptor_filename = NULL;
    const char *ligand_filename = NULL;
    if (!PyArg_ParseTuple(args, "|ss", &receptor_filename, &ligand_filename)) {
        Err("Parameter type error!");
        return 0;
    }

    if (receptor_filename == NULL || ligand_filename == NULL) {
        return 0;
    }

    // grid
    PyObject *box_object = PyDict_GetItemString(kwargs, "box");
    double box = (box_object == NULL) ? 10.0 : o2d(box_object);
    
    PyObject *bin_object = PyDict_GetItemString(kwargs, "bin");
    double bin = (bin_object == NULL) ? 0.2 : o2d(bin_object);

    auto &&rec = jnc::pdb::read_pdb(receptor_filename);
    auto &&lig = jnc::mol2::read_mol2s(ligand_filename)[0];
    auto &&center = residue_center(lig.atoms);

    PocketGrid *grid = new PocketGrid;
    grid->bin = bin;
    grid->box = box;
    grid->center = center;
    grid->size = int(box / bin);

    self->pocket->grid = grid;

    // receptor
    std::list<Atom> atoms;
    for (auto && chain : rec[0]) {
        for (auto && res : chain) {
            for (auto && atom : res) {
                if (atom_in_box(atom, center, box+6.0)) {
                    atoms.push_back(atom);
                }
            }
        }
    }

    PocketReceptor *p = new PocketReceptor;

    int natoms = atoms.size();
    p->atoms.resize(natoms);
    int iatom = 0;
    for (auto && atom : atoms) {
        auto &pa = p->atoms[iatom];
        for (int i = 0; i < 3; i++) {
            pa.coord[i] = atom[i];
            pa.index[i] = int((atom[i] - (center[i] - box/2.0)) / bin);
        }
        pa.isLigand = 0;
        pa.name = atom.name;
        iatom++;
    }

    self->pocket->receptor = p;

    // id
    self->pocket->id = 0;
    jnc::hash_combine(self->pocket->id, std::string(receptor_filename));
    jnc::hash_combine(self->pocket->id, std::string(ligand_filename));

    return 0;
    // return Py_BuildValue("{s:i,s:N,s:d,s:N}", "id", int(id), "center", a2o(center, 3), "box", box, "atoms", atoms_object);//     return 0;
}

static void Pocket_dealloc(PocketObject *self) {
    // std::cout << "dealloc: " << self->pocket->id << std::endl;
    if (self->pocket != NULL) {
        if (self->pocket->ligand != NULL) {
            delete self->pocket->ligand;
        } else {
            if (self->pocket->receptor != NULL) {
                delete self->pocket->grid;
                delete self->pocket->receptor;
            }
        }
        delete self->pocket;
    }
    Py_TYPE(self)->tp_free((PyObject *) self);
}

std::vector<std::array<double, 3>> fibonacci_sphere(double radius = 1, int samples=100) {
    std::vector<std::array<double, 3>> points(samples);
    double phi = jnc::pi * (3. - std::sqrt(5.));  // golden angle in radians

    for (int i = 0; i < samples; i++) {
        double y = (1 - (i / float(samples - 1)) * 2) * radius;  // y goes from 1 to -1
        double r = std::sqrt(radius * radius - y * y);  // radius at y

        double theta = phi * i;  // golden angle increment

        double x = std::cos(theta) * r;
        double z = std::sin(theta) * r;

        points[i][0] = x;
        points[i][1] = y;
        points[i][2] = z;
    }

    return points;
}

static std::vector<std::array<double, 3>> bond_points = fibonacci_sphere(1, 50);

template<typename Index_>
static std::list<std::array<int, 3>> surrounding_positions(const PocketGrid &grid, const Index_ &ind, double r1, double r2) {
    std::list<std::array<int, 3>> positions;

    double halfBin = grid.bin / 2.0;
    double halfBox = grid.box / 2.0;

    double ox = grid.center[0] - halfBox;
    double oy = grid.center[1] - halfBox;
    double oz = grid.center[2] - halfBox;
    // std::cout << "bin: " << grid.bin << ", box: " << grid.box << ", [xyz]: " << ox << ' ' << oy << ' ' << oz << std::endl;

    double x = grid.bin * ind[0] + ox;
    double y = grid.bin * ind[1] + oy;
    double z = grid.bin * ind[2] + oz;

    int x1 = (int) ((x - r2 - ox) / grid.bin);
    int y1 = (int) ((y - r2 - oy) / grid.bin);
    int z1 = (int) ((z - r2 - oz) / grid.bin);

    int x2 = (int) ((x + r2 - ox) / grid.bin);
    int y2 = (int) ((y + r2 - oy) / grid.bin);
    int z2 = (int) ((z + r2 - oz) / grid.bin);

    // std::cout << "bin: " << grid.bin << ", box: " << grid.box << ", [xyz]: " << ox << ' ' << oy << ' ' << oz << std::endl;
    // std::cout << x << ' ' << y << ' ' << z << ' ' << x1 << ' ' << y1 << ' ' << z1 << ' ' << x2 << ' ' << y2 << ' ' << z2 << std::endl;

    for (int i = x1; i <= x2; i++) {
        if (i < 0 || i > grid.size) continue;
        for (int j = y1; j <= y2; j++) {
            if (j < 0 || j > grid.size) continue;
            for (int k = z1; k <= z2; k++) {
                if (k < 0 || k > grid.size) continue;
                double dx = i * grid.bin + ox - x;
                double dy = j * grid.bin + oy - y;
                double dz = k * grid.bin + oz - z;
                double d2 = dx * dx + dy * dy + dz * dz;

                double l = r1 - halfBin;
                if (l < 0) l = 0;

                double u = r2 + halfBin;
                
                if (d2 >= l * l && d2 < u * u) {
                    positions.push_back({i, j, k});
                }
            }
        }
    }
    return positions;
}

static PocketObject *make_pocket();

static PyObject *Pocket_find_children(PocketObject *self, PyObject *args, PyObject *kwargs) {
    std::list<std::array<int, 3>> positions;
    PocketLigandAtom *prevAtom;
    PocketGrid *grid = self->pocket->grid;
    std::unordered_set<std::array<int, 3>, jnc::ArrayHash<int, 3>> occupied;
    std::unordered_set<std::array<int, 3>, jnc::ArrayHash<int, 3>> existed;
    for (auto && atom : self->pocket->receptor->atoms) {
        for (auto && ind : surrounding_positions(*grid, atom.index, 0, 2)) { // hydrogen bond length: 2.5-4
            occupied.insert(ind);
        }
    }

    if (self->pocket->ligand == NULL) { // root
        prevAtom = NULL;
        for (auto && atom : self->pocket->receptor->atoms) {
            for (auto && ind : surrounding_positions(*grid, atom.index, 2, 4)) { // hydrogen bond length: 2.5-4
                if (existed.find(ind) == existed.end()) {
                    positions.push_back(ind);
                    existed.insert(ind);
                }
            }
        }
    } else { // non-root
        prevAtom = &(self->pocket->ligand->lastAtom);
        auto *atom = prevAtom;
        do {
            for (auto && ind : surrounding_positions(*grid, atom->index, 0, 1)) { // bond length: 1-1.5
                occupied.insert(ind);
            }
            for (auto && ind : surrounding_positions(*grid, atom->index, 1, 1.5)) { // bond length: 1-1.5
                if (existed.find(ind) == existed.end()) {
                    existed.insert(ind);
                    positions.push_back(ind);
                }
            }
            atom = atom->prev;
        } while (atom != NULL);
    }
    // std::cout << positions.size() << " positions" << std::endl;

    PyObject *children = PyList_New(0);
    for (auto && ind : positions) {
        if (occupied.find(ind) == occupied.end()) {
            auto child = make_pocket();
            // std::cout << "Child reference count1: " << Py_REFCNT(child) << std::endl;
            child->pocket->receptor = self->pocket->receptor;
            child->pocket->grid = self->pocket->grid;
            child->pocket->id = self->pocket->id;
            for (int i = 0; i < 3; i++) jnc::hash_combine(child->pocket->id, ind[i]);

            child->pocket->ligand = new PocketLigand;

            auto &pa = child->pocket->ligand->lastAtom;
            pa.element = 0;
            for (int i = 0; i < 3; i++) pa.index[i] = ind[i];
            pa.prev = prevAtom;

            PyList_Append(children, (PyObject *) child);
            // std::cout << "Child reference count2: " << Py_REFCNT(child) << std::endl;
            Py_DECREF(child);
            // std::cout << "Child reference count3: " << Py_REFCNT(child) << std::endl;
        }
    }

    return children;
}

static PyObject *Pocket_receptor_grid(PocketObject *self, PyObject *args, PyObject *kwargs) {
    auto &atoms = self->pocket->receptor->atoms;
    PyObject *ls = PyList_New(0);
    auto *grid = self->pocket->grid;
    for (auto && atom : atoms) {
        if (std::all_of(atom.index.begin(), atom.index.end(), [grid](int i){
            return i >= 0 && i <= grid->size;
        })) {
            PyObject *ind = PyList_New(3);
            for (int i = 0; i < 3; i++) PyList_SET_ITEM(ind, i, PyLong_FromLong(atom.index[i]));
            PyList_Append(ls, ind);
            Py_DECREF(ind);
        }
    }
    return ls;
}

static PyObject *Pocket_ligand_grid(PocketObject *self, PyObject *args, PyObject *kwargs) {
    PyObject *ls = PyList_New(0);
    if (self->pocket->ligand != NULL) {
        auto atom = &(self->pocket->ligand->lastAtom);
        do {
            PyObject *ind = PyList_New(3);
            for (int i = 0; i < 3; i++) PyList_SET_ITEM(ind, i, PyLong_FromLong(atom->index[i]));
            PyList_Append(ls, ind);
            Py_DECREF(ind);
            atom = atom->prev;
        } while (atom != NULL);
    }
    return ls;
}

static PyMethodDef Pocket_methods[] = {
    {"find_children", (PyCFunction) Pocket_find_children, METH_NOARGS, "find children"},
    {"receptor_grid", (PyCFunction) Pocket_receptor_grid, METH_NOARGS, "receptor grid"},
    {"ligand_grid", (PyCFunction) Pocket_ligand_grid, METH_NOARGS, "ligand grid"},
    {NULL}  /* Sentinel */
};

static PyTypeObject PocketType = []{
    PyTypeObject obj{PyVarObject_HEAD_INIT(NULL, 0)};
    obj.tp_name = "Pocket",
    obj.tp_basicsize = sizeof(PocketObject),
    obj.tp_itemsize = 0,
    obj.tp_dealloc = (destructor) Pocket_dealloc,
    obj.tp_flags = Py_TPFLAGS_DEFAULT;
    obj.tp_doc = "Pocket Object";
    obj.tp_init = (initproc) Pocket_init;
    obj.tp_new = Pocket_new;
    obj.tp_methods = Pocket_methods;
    return obj;
}();

static PocketObject *make_pocket() {
    PocketObject *obj = (PocketObject *) PyObject_CallObject((PyObject *) &PocketType, NULL);
    return obj;
}

//class Mat3 {
//public:
//    double *data;
//    std::array<int, 3> shape;
//    std::array<int, 3> strides;
//
//    const Mat3 &fft(Mat3 &&mat) const {
//        fftw_plan plan = fftw_plan_dft_3d(shape[0], shape[1], shape[2], reinterpret_cast<fftw_complex*>(data), reinterpret_cast<fftw_complex*>(mat.data), FFTW_FORWARD, FFTW_ESTIMATE);
//        fftw_execute(plan); /* repeat as needed */
//        fftw_destroy_plan(plan);
//        return *this;
//    }
//
//    const Mat3 &ifft(Mat3 &&mat) const {
//        fftw_plan plan = fftw_plan_dft_3d(shape[0], shape[1], shape[2], reinterpret_cast<fftw_complex*>(data), reinterpret_cast<fftw_complex*>(mat.data), FFTW_BACKWARD, FFTW_ESTIMATE);
//        fftw_execute(plan); /* repeat as needed */
//        fftw_destroy_plan(plan);
//        return *this;
//    }
//};

//static PyObject *amorphize(PyObject *atoms, int n) {
//    Mat3 radius_grid, radius_fgrid, grid1, fgrid1, grid2, fgrid2;
//
//    PyObject *ls = PyList_New(0);
//    for (int i = 0; i < n; i++) {
//        auto as = copy_atoms(atoms);
//        set_grid(grid2, atoms);
//        score = fit(grid1, grid2);
//        update_grid(grid2, as);
//        PyList_Append(ls, as);
//    }
//    return ls
//}
//
//static PyObject *amorphize(PyObject *self, PyObject *args) {
//    PyObject *atoms;
//    int n;
//    PARSE_TUPLE(args, "Oi", &atoms, &n);
//    return amorphize_(atoms, n);
//}

static PyMethodDef jnpy_methods[] = {
    {"add", add, METH_VARARGS, "Add"},
//    {"newatom", newatom, METH_VARARGS, "New Atom"},
    {"suppos", suppos, METH_VARARGS, "Superposion"},
    {"sp_apply", sp_apply, METH_VARARGS, "Apply Superposion"},
    {"align", align, METH_VARARGS, "Apply Superposion"},
    {"rmsd", rmsd, METH_VARARGS, "RMSD"},
    {"map_atoms", map_atoms, METH_VARARGS, "Map Atoms"},
//    {"Amorphize", amorphize, METH_VARARGS, "Amorphize"},
//    {"aa321", aa321, METH_VARARGS, "aa321"},
//    {"aa123", aa123, METH_VARARGS, "aa123"},
    {NULL, NULL, 0, NULL} // sentinel
};

// static struct PyModuleDef jnpy_module = {
//     PyModuleDef_HEAD_INIT,
//     "_jnpy",   /* name of module */
//     NULL, /* module documentation, may be NULL */
//     -1,       /* size of per-interpreter state of the module,
//                  or -1 if the module keeps state in global variables. */
//     jnpy_methods
// };

static PyModuleDef jnpy_module = {
    PyModuleDef_HEAD_INIT,
    .m_name = "_jnpy",
    .m_doc = "Jian's Python Library",
    .m_size = -1,
    .m_methods = jnpy_methods,
};

PyMODINIT_FUNC PyInit__jnpy(void) {
    // return PyModule_Create(&jnpy_module);

    PyObject *m;

    if (PyType_Ready(&PocketType) < 0)
        return NULL;

    m = PyModule_Create(&jnpy_module);
    if (m == NULL)
        return NULL;

    Py_INCREF(&PocketType);
    if (PyModule_AddObject(m, "Pocket", (PyObject *) &PocketType) < 0) {
        Py_DECREF(&PocketType);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}

/*
PyObject *function = PyObject_GetAttrString(add_module, "add");
PyObject *args = PyTuple_New(0);
PyObject *kwargs = Py_BuildValue("{s:i}", "b", 5)
result = PyObject_Call(function, args, kwargs);

Py_DECREF(kwargs);
Py_DECREF(args);
Py_DECREF(function);
*/

/*
PyObject *module = PyImport_AddModule( "gto" );
PyObject *moduleDict = PyModule_GetDict( module );

PyObject *classObj = PyDict_GetItemString( moduleDict, "ObjectInfo" );

PyObject *args = Py_BuildValue( "()" ); // Empty tuple
PyObject *objInfo = PyInstance_New( classObj, args, NULL );

PyObject_SetAttr( objInfo, 
PyString_FromString( "name" ),
PyString_FromString( 
reader->stringFromId( oi.name ).c_str() ) );
*/

/*
 * Build a dictionary
 */
// Py_BuildValue("{s:i,s:O}", "abc", a_c_int, "def", a_python_list);  


