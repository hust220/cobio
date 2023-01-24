#include "./py_pdb.h"
#include "./jnc_geom.h"
#include "./jnc_bio_pdb.h"

using Residue = jnc::bio::PdbResidue;
using Atom = jnc::bio::PdbAtom;

Residue obj2residue(PyObject *o) {
  Residue residue;
  return std::move(residue);
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

typedef struct {
  PyObject_HEAD jnc::bio::Pdb *pdb;
} PdbObject;

PyObject *Pdb_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  PdbObject *self;
  self = (PdbObject *)type->tp_alloc(type, 0);
  if (self != NULL) {
    self->pdb = NULL;
  }
  return (PyObject *)self;
}

int Pdb_init(PdbObject *self, PyObject *args, PyObject *kwargs) {
  // const char *filename = NULL;
  // if (!PyArg_ParseTuple(args, "s", &filename)) {
  //   Err("Parameter type error!");
  //   return 1;
  // }
  int sz = PyTuple_GET_SIZE(args);
  auto obj = PyTuple_GET_ITEM(args, 0);

  if (PyUnicode_Check(obj)) {
    self->pdb = new jnc::bio::Pdb(o2s(obj));
    std::cout << o2s(obj) << std::endl;
  } else {
    self->pdb = new jnc::bio::Pdb;
    jnc::bio::PdbReader reader(*self->pdb);
    while (true) {
      auto r = PyFile_GetLine(obj, 0);
      auto &&line = std::string(o2s(r));
      if (line.empty()) {
        Py_DECREF(r);
        break;
      } else {
        reader.read_line(line);
        Py_DECREF(r);
      }
    }
    reader.end_reading();
  }

  //   std::cout << self->pdb->size() << std::endl;

  return 0;
}

void Pdb_dealloc(PdbObject *self) {
  if (self->pdb != NULL) {
    delete self->pdb;
  }
  Py_TYPE(self)->tp_free((PyObject *)self);
}

PyObject *Pdb_seq(PdbObject *self, PyObject *args, PyObject *kwargs) {
  std::string s;
  auto &m = jnc::bio::g_aa3;
  for (auto &&chain : self->pdb->at(0)) {
    for (auto &&res : chain) {
      if (m.find(res.name) != m.end()) {
        s += m.at(res.name);
      } else {
        s += 'X';
      }
    }
  }
  return str2obj(s.c_str());
}

PyObject *Pdb_contacts(PdbObject *self, PyObject *args, PyObject *kwargs) {
  PyObject *cutoff_object = NULL;
  if (kwargs != NULL) {
    cutoff_object = PyDict_GetItemString(kwargs, "cutoff");
  }
  double cutoff = (cutoff_object == NULL) ? 10.0 : o2d(cutoff_object);
  double cutoff2 = cutoff * cutoff;
  double bin = cutoff;

  auto &&rs = self->pdb->at(0).presidues();
  int nres = rs.size();
  std::vector<std::array<double, 3>> coords(nres);
  std::array<double, 3> origin = {0, 0, 0};
  for (int i = 0; i < nres; i++) {
    auto &r = *rs[i];
    auto &c = coords[i];
    if (r.has_atom("CA")) {
      auto &atom = r.atom("CA");
      for (int j = 0; j < 3; j++)
        c[j] = atom[j];
    } else {
      for (int j = 0; j < 3; j++)
        c[j] = 0;
      for (auto &&atom : r) {
        for (int j = 0; j < 3; j++)
          c[j] += atom[j];
      }
      for (int j = 0; j < 3; j++)
        c[j] /= r.size();
    }
    for (int j = 0; j < 3; j++) {
      if (i == 0 || origin[j] > c[j]) {
        origin[j] = c[j];
      }
    }
  }

  std::unordered_map<std::array<int, 3>, std::list<int>, jnc::ArrayHash<int, 3>> grid;
  std::vector<std::array<int, 3>> indices(nres);
  for (int i = 0; i < nres; i++) {
    for (int j = 0; j < 3; j++) {
      indices[i][j] = int((coords[i][j] - origin[j]) / bin);
    }
    grid[indices[i]].push_back(i);
  }

  std::list<std::tuple<int, int, double>> contacts;
  for (int i = 0; i < nres; i++) {
    auto &ind = indices[i];

    for (int a = ind[0] - 1; a <= ind[0] + 1; a++) {
      for (int b = ind[1] - 1; b <= ind[1] + 1; b++) {
        for (int c = ind[2] - 1; c <= ind[2] + 1; c++) {
          auto &neighbors = grid[std::array<int, 3>{a, b, c}];
          for (auto &j : neighbors) {
            if (i < j) {
              double d2 = jnc::geom::dist2(coords[i], coords[j]);
              if (d2 < cutoff2) {
                contacts.push_back(std::make_tuple(i, j, std::sqrt(d2)));
              }
            }
          }
        }
      }
    }
  }

  int ncontacts = contacts.size();
  PyObject *contacts_object = PyList_New(ncontacts);
  int icontact = 0;
  int ir, jr;
  double d;
  for (auto &contact : contacts) {
    std::tie(ir, jr, d) = contact;
    PyObject *contact_object = PyList_New(3);
    PyList_SET_ITEM(contact_object, 0, PyLong_FromLong(ir));
    PyList_SET_ITEM(contact_object, 1, PyLong_FromLong(jr));
    PyList_SET_ITEM(contact_object, 2, PyFloat_FromDouble(d));
    PyList_SET_ITEM(contacts_object, icontact, contact_object);
    icontact++;
  }

  return contacts_object;
}

PyMethodDef Pdb_methods[] = {
    {"seq", (PyCFunction)Pdb_seq, METH_NOARGS, "Get the sequence from a protein file"},
    {"contacts", (PyCFunction)Pdb_contacts, METH_VARARGS | METH_KEYWORDS, "Get the contacts from a protein file"},
    {NULL} /* Sentinel */
};

PyTypeObject PdbType = [] {
  PyTypeObject obj{PyVarObject_HEAD_INIT(NULL, 0)};
  obj.tp_name = "Pdb", obj.tp_basicsize = sizeof(PdbObject), obj.tp_itemsize = 0,
  obj.tp_dealloc = (destructor)Pdb_dealloc, obj.tp_flags = Py_TPFLAGS_DEFAULT;
  obj.tp_doc = "Pdb Object";
  obj.tp_init = (initproc)Pdb_init;
  obj.tp_new = Pdb_new;
  obj.tp_methods = Pdb_methods;
  return obj;
}();

