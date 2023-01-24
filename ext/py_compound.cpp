#include "./py_compound.h"
#include "./jnc_bio_pdb.h"
#include "./jnc_bio_mol2.h"

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

struct PocketGrid {
  std::array<double, 3> center;
  double box;
  double bin;
  int size;
};

template <typename Grid_, typename Coord_> std::array<int, 3> grid_c2i(const Grid_ &grid, const Coord_ &c) {
  std::array<int, 3> ind;
  for (int i = 0; i < 3; i++)
    ind[i] = int((c[i] - (grid.center[i] - grid.box / 2.0)) / grid.bin);
  return ind;
}

template <typename Grid_, typename Ind_> std::array<double, 3> grid_i2c(const Grid_ &grid, const Ind_ &ind) {
  std::array<double, 3> c;
  for (int i = 0; i < 3; i++)
    c[i] = ind[i] * grid.bin + (grid.center[i] - grid.box / 2.0);
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
  PocketGrid *grid = NULL;         // shared by all the nodes
  PocketReceptor *receptor = NULL; // shared by all the nodes
  PocketLigand *ligand = NULL;     // owned by the current node
};

template <typename Atom1_, typename Atom2_> bool atom_in_box(Atom1_ &&atom1, Atom2_ &&atom2, double box) {
  double halfBox = box / 2.0;
  double x = atom1[0] - atom2[0];
  double y = atom1[1] - atom2[1];
  double z = atom1[2] - atom2[2];
  return std::abs(x) < halfBox && std::abs(y) < halfBox && std::abs(z) < halfBox;
}

template <typename Residue_> std::array<double, 3> residue_center(const Residue_ &residue) {
  std::array<double, 3> c{0, 0, 0};
  int iatom = 0;
  for (auto &&atom : residue) {
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
  PyObject_HEAD Pocket *pocket;
} PocketObject;

PocketObject *make_pocket();

// static std::vector<std::array<double, 3>> bond_points = fibonacci_sphere(1, 50);

template <typename Index_>
static std::list<std::array<int, 3>> surrounding_positions(const PocketGrid &grid, const Index_ &ind, double r1,
                                                           double r2) {
  std::list<std::array<int, 3>> positions;

  double halfBin = grid.bin / 2.0;
  double halfBox = grid.box / 2.0;

  double ox = grid.center[0] - halfBox;
  double oy = grid.center[1] - halfBox;
  double oz = grid.center[2] - halfBox;
  // std::cout << "bin: " << grid.bin << ", box: " << grid.box << ", [xyz]: " <<
  // ox << ' ' << oy << ' ' << oz << std::endl;

  double x = grid.bin * ind[0] + ox;
  double y = grid.bin * ind[1] + oy;
  double z = grid.bin * ind[2] + oz;

  int x1 = (int)((x - r2 - ox) / grid.bin);
  int y1 = (int)((y - r2 - oy) / grid.bin);
  int z1 = (int)((z - r2 - oz) / grid.bin);

  int x2 = (int)((x + r2 - ox) / grid.bin);
  int y2 = (int)((y + r2 - oy) / grid.bin);
  int z2 = (int)((z + r2 - oz) / grid.bin);

  // std::cout << "bin: " << grid.bin << ", box: " << grid.box << ", [xyz]: " <<
  // ox << ' ' << oy << ' ' << oz << std::endl; std::cout << x << ' ' << y << '
  // ' << z << ' ' << x1 << ' ' << y1 << ' ' << z1 << ' ' << x2 << ' ' << y2 <<
  // ' ' << z2 << std::endl;

  for (int i = x1; i <= x2; i++) {
    if (i < 0 || i > grid.size)
      continue;
    for (int j = y1; j <= y2; j++) {
      if (j < 0 || j > grid.size)
        continue;
      for (int k = z1; k <= z2; k++) {
        if (k < 0 || k > grid.size)
          continue;
        double dx = i * grid.bin + ox - x;
        double dy = j * grid.bin + oy - y;
        double dz = k * grid.bin + oz - z;
        double d2 = dx * dx + dy * dy + dz * dz;

        double l = r1 - halfBin;
        if (l < 0)
          l = 0;

        double u = r2 + halfBin;

        if (d2 >= l * l && d2 < u * u) {
          positions.push_back({i, j, k});
        }
      }
    }
  }
  return positions;
}

PyObject *Pocket_find_children(PocketObject *self, PyObject *args, PyObject *kwargs);
PyObject *Pocket_receptor_grid(PocketObject *self, PyObject *args, PyObject *kwargs);
PyObject *Pocket_ligand_grid(PocketObject *self, PyObject *args, PyObject *kwargs);


PyObject *compound_grow(PyObject *self, PyObject *args, PyObject *kwargs) {
  PyObject *atoms;
  double bin = 0.4;
  int ntypes = 7;
  int coords_only = 0;

  char *kwlist[] = {"atoms", "bin", "ntypes", "coords_only", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|dip", kwlist, &atoms, &bin, &ntypes, &coords_only)) {
    Err("Parameter type error!");
    return NULL;
  }

  int natoms = PyObject_Size(atoms);

  if (natoms == 0) {
    if (coords_only) {
      PyObject *children_object = PyList_New(1);
      PyList_SetItem(children_object, 0, Py_BuildValue("(i,i,i)", 10, 10, 10)); // steal the reference
      return children_object;
    } else {
      PyObject *children_object = PyList_New(ntypes);
      for (int i = 0; i < ntypes; i++) {
        PyList_SetItem(children_object, i, Py_BuildValue("(i,i,i,i)", 10, 10, 10, i)); // steal the reference
      }
      return children_object;
    }
  }

  std::unordered_set<std::array<int, 3>, jnc::ArrayHash<int, 3>> occupied;
  double min_bond = 0.8;
  double max_bond = 1.2;
  int t = int(min_bond / bin - 1);
  int r = int(std::ceil(max_bond / bin - 1));
  std::vector<std::array<int, 3>> coords(natoms);
  for (int iatom = 0; iatom < natoms; iatom++) {
    auto atom = obj_get<PyObject *>(atoms, iatom); // borrowed reference

    auto &c = coords[iatom];
    for (int i = 0; i < 3; i++) {
      c[i] = obj_get<int>(atom, i);
    }
    for (int i = c[0] - t; i <= c[0] + t; i++) {
      for (int j = c[1] - t; j <= c[1] + t; j++) {
        for (int k = c[2] - t; k <= c[2] + t; k++) {
          double x = (std::abs(i - c[0]) + 1) * bin;
          double y = (std::abs(j - c[1]) + 1) * bin;
          double z = (std::abs(k - c[2]) + 1) * bin;
          if (x * x + y * y + z * z < min_bond * min_bond) {
            occupied.insert(std::array<int, 3>{i, j, k});
          }
        }
      }
    }
    // int type = PyLong_AS_LONG(PyList_GET_ITEM(atom, 3));
  }

  std::unordered_set<std::array<int, 3>, jnc::ArrayHash<int, 3>> children;
  for (int iatom = 0; iatom < natoms; iatom++) {
    auto &c = coords[iatom];
    for (int i = c[0] - r; i <= c[0] + r; i++) {
      for (int j = c[1] - r; j <= c[1] + r; j++) {
        for (int k = c[2] - r; k <= c[2] + r; k++) {
          std::array<int, 3> ind{i, j, k};
          double x = std::max(std::abs(i - c[0]) - 1, 0) * bin;
          double y = std::max(std::abs(j - c[1]) - 1, 0) * bin;
          double z = std::max(std::abs(k - c[2]) - 1, 0) * bin;
          if (x * x + y * y + z * z < max_bond * max_bond) {
            if (occupied.find(ind) == occupied.end()) {
              if (i >= 0 && j >= 0 && k >= 0) {
                children.insert(ind);
              }
            }
          }
        }
      }
    }
  }

  int nchildren = children.size() * (coords_only ? 1 : ntypes);

  PyObject *children_object = PyList_New(nchildren);
  int ichild = 0;
  for (auto &&c : children) {
    if (coords_only) {
      PyList_SetItem(children_object, ichild, Py_BuildValue("(i,i,i)", c[0], c[1], c[2])); // steal the reference
      ichild++;
    } else {
      for (int i = 0; i < ntypes; i++) {
        PyList_SetItem(children_object, ichild, Py_BuildValue("(i,i,i,i)", c[0], c[1], c[2], i)); // steal the reference
        ichild++;
      }
    }
  }

  return children_object;
}

PyObject *Pocket_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  PocketObject *self;
  self = (PocketObject *)type->tp_alloc(type, 0);
  if (self != NULL) {
    self->pocket = new Pocket;
    self->pocket->receptor = NULL;
    self->pocket->receptor = NULL;
    self->pocket->ligand = NULL;
    self->pocket->id = 0;
  }
  return (PyObject *)self;
}

int Pocket_init(PocketObject *self, PyObject *args, PyObject *kwargs) {
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

  auto &&rec = jnc::bio::read_pdb(receptor_filename);
  auto &&lig = jnc::bio::read_mol2s(ligand_filename)[0];
  auto &&center = residue_center(lig.atoms);

  PocketGrid *grid = new PocketGrid;
  grid->bin = bin;
  grid->box = box;
  grid->center = center;
  grid->size = int(box / bin);

  self->pocket->grid = grid;

  // receptor
  std::list<jnc::bio::PdbAtom> atoms;
  for (auto &&chain : rec[0]) {
    for (auto &&res : chain) {
      for (auto &&atom : res) {
        if (atom_in_box(atom, center, box + 6.0)) {
          atoms.push_back(atom);
        }
      }
    }
  }

  PocketReceptor *p = new PocketReceptor;

  int natoms = atoms.size();
  p->atoms.resize(natoms);
  int iatom = 0;
  for (auto &&atom : atoms) {
    auto &pa = p->atoms[iatom];
    for (int i = 0; i < 3; i++) {
      pa.coord[i] = atom[i];
      pa.index[i] = int((atom[i] - (center[i] - box / 2.0)) / bin);
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
  // return Py_BuildValue("{s:i,s:N,s:d,s:N}", "id", int(id), "center",
  // a2o(center, 3), "box", box, "atoms", atoms_object);//     return 0;
}

void Pocket_dealloc(PocketObject *self) {
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
  Py_TYPE(self)->tp_free((PyObject *)self);
}

std::vector<std::array<double, 3>> fibonacci_sphere(double radius = 1, int samples = 100) {
  std::vector<std::array<double, 3>> points(samples);
  double phi = jnc::pi * (3. - std::sqrt(5.)); // golden angle in radians

  for (int i = 0; i < samples; i++) {
    double y = (1 - (i / float(samples - 1)) * 2) * radius; // y goes from 1 to -1
    double r = std::sqrt(radius * radius - y * y);          // radius at y

    double theta = phi * i; // golden angle increment

    double x = std::cos(theta) * r;
    double z = std::sin(theta) * r;

    points[i][0] = x;
    points[i][1] = y;
    points[i][2] = z;
  }

  return points;
}

PyObject *Pocket_find_children(PocketObject *self, PyObject *args, PyObject *kwargs) {
  std::list<std::array<int, 3>> positions;
  PocketLigandAtom *prevAtom;
  PocketGrid *grid = self->pocket->grid;
  std::unordered_set<std::array<int, 3>, jnc::ArrayHash<int, 3>> occupied;
  std::unordered_set<std::array<int, 3>, jnc::ArrayHash<int, 3>> existed;
  for (auto &&atom : self->pocket->receptor->atoms) {
    for (auto &&ind : surrounding_positions(*grid, atom.index, 0,
                                            2)) { // hydrogen bond length: 2.5-4
      occupied.insert(ind);
    }
  }

  if (self->pocket->ligand == NULL) { // root
    prevAtom = NULL;
    for (auto &&atom : self->pocket->receptor->atoms) {
      for (auto &&ind : surrounding_positions(*grid, atom.index, 2, 4)) { // hydrogen bond length: 2.5-4
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
      for (auto &&ind : surrounding_positions(*grid, atom->index, 0,
                                              1)) { // bond length: 1-1.5
        occupied.insert(ind);
      }
      for (auto &&ind : surrounding_positions(*grid, atom->index, 1,
                                              1.5)) { // bond length: 1-1.5
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
  for (auto &&ind : positions) {
    if (occupied.find(ind) == occupied.end()) {
      auto child = make_pocket();
      // std::cout << "Child reference count1: " << Py_REFCNT(child) <<
      // std::endl;
      child->pocket->receptor = self->pocket->receptor;
      child->pocket->grid = self->pocket->grid;
      child->pocket->id = self->pocket->id;
      for (int i = 0; i < 3; i++)
        jnc::hash_combine(child->pocket->id, ind[i]);

      child->pocket->ligand = new PocketLigand;

      auto &pa = child->pocket->ligand->lastAtom;
      pa.element = 0;
      for (int i = 0; i < 3; i++)
        pa.index[i] = ind[i];
      pa.prev = prevAtom;

      PyList_Append(children, (PyObject *)child);
      // std::cout << "Child reference count2: " << Py_REFCNT(child) <<
      // std::endl;
      Py_DECREF(child);
      // std::cout << "Child reference count3: " << Py_REFCNT(child) <<
      // std::endl;
    }
  }

  return children;
}

PyObject *Pocket_receptor_grid(PocketObject *self, PyObject *args, PyObject *kwargs) {
  auto &atoms = self->pocket->receptor->atoms;
  PyObject *ls = PyList_New(0);
  auto *grid = self->pocket->grid;
  for (auto &&atom : atoms) {
    if (std::all_of(atom.index.begin(), atom.index.end(), [grid](int i) { return i >= 0 && i <= grid->size; })) {
      PyObject *ind = PyList_New(3);
      for (int i = 0; i < 3; i++)
        PyList_SET_ITEM(ind, i, PyLong_FromLong(atom.index[i]));
      PyList_Append(ls, ind);
      Py_DECREF(ind);
    }
  }
  return ls;
}

PyObject *Pocket_ligand_grid(PocketObject *self, PyObject *args, PyObject *kwargs) {
  PyObject *ls = PyList_New(0);
  if (self->pocket->ligand != NULL) {
    auto atom = &(self->pocket->ligand->lastAtom);
    do {
      PyObject *ind = PyList_New(3);
      for (int i = 0; i < 3; i++)
        PyList_SET_ITEM(ind, i, PyLong_FromLong(atom->index[i]));
      PyList_Append(ls, ind);
      Py_DECREF(ind);
      atom = atom->prev;
    } while (atom != NULL);
  }
  return ls;
}

PyMethodDef Pocket_methods[] = {
    {"find_children", (PyCFunction)Pocket_find_children, METH_NOARGS, "find children"},
    {"receptor_grid", (PyCFunction)Pocket_receptor_grid, METH_NOARGS, "receptor grid"},
    {"ligand_grid", (PyCFunction)Pocket_ligand_grid, METH_NOARGS, "ligand grid"},
    {NULL} /* Sentinel */
};

PyTypeObject PocketType = [] {
  PyTypeObject obj{PyVarObject_HEAD_INIT(NULL, 0)};
  obj.tp_name = "Pocket", obj.tp_basicsize = sizeof(PocketObject), obj.tp_itemsize = 0,
  obj.tp_dealloc = (destructor)Pocket_dealloc, obj.tp_flags = Py_TPFLAGS_DEFAULT;
  obj.tp_doc = "Pocket Object";
  obj.tp_init = (initproc)Pocket_init;
  obj.tp_new = Pocket_new;
  obj.tp_methods = Pocket_methods;
  return obj;
}();

PocketObject *make_pocket() {
  PocketObject *obj = (PocketObject *)PyObject_CallObject((PyObject *)&PocketType, NULL);
  return obj;
}

