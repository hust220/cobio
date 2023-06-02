#include "py_pdb.h"
#include "./jnc_bio_pdb.h"
#include "./jnc_geom.h"
#include "./py_cif.h"

using Residue = jnc::bio::PdbResidue;
using Atom = jnc::bio::PdbAtom;

struct cif_line_t {
  int atom_num;// id: 
  std::string atom_element;// type_symbol: 
  std::string atom_name;// label_atom_id: 
  std::string residue_name;// label_comp_id: 
  std::string chain_name;// label_asym_id: 
  int residue_num;// label_seq_id: 
  double x;// Cartn_x: 
  double y;// Cartn_y: 
  double z;// Cartn_z: 
  int model_num; // pdbx_PDB_model_num
};

void parse_cif_line(cif_line_t &line, const std::map<std::string, std::string> &values) {
  for (auto &&p : values) {
    if (p.first == "id") line.atom_num = JN_INT(p.second);
    else if (p.first == "label_atom_id") line.atom_name = p.second;
    else if (p.first == "label_comp_id") line.residue_name = p.second;
    else if (p.first == "label_asym_id") line.chain_name = p.second;
    else if (p.first == "label_seq_id") line.residue_num = JN_INT(p.second);
    else if (p.first == "Cartn_x") line.x = JN_DBL(p.second);
    else if (p.first == "Cartn_y") line.y = JN_DBL(p.second);
    else if (p.first == "Cartn_z") line.z = JN_DBL(p.second);
    else if (p.first == "pdbx_PDB_model_num") line.model_num = JN_INT(p.second);
  }
}

PyObject *new_atom() {
  PyObject *dict = PyDict_New();
  PyDict_SetItemString(dict, "name", Py_BuildValue("s", ""));
  PyDict_SetItemString(dict, "num", Py_BuildValue("i", 0));
  PyDict_SetItemString(dict, "coords", PyList_New(3));
  PyDict_SetItemString(dict, "is_het", PyBool_FromLong(0));
  return dict;
}

PyObject *new_residue() {
  PyObject *dict = PyDict_New();
  PyDict_SetItemString(dict, "num", Py_BuildValue("i", 0));
  PyDict_SetItemString(dict, "name", Py_BuildValue("s", ""));
  PyDict_SetItemString(dict, "atoms", PyList_New(0));
  return dict;
}

PyObject *new_chain() {
  PyObject *dict = PyDict_New();
  PyDict_SetItemString(dict, "name", Py_BuildValue("s", ""));
  PyDict_SetItemString(dict, "residues", PyList_New(0));
  return dict;
}

PyObject *new_model() {
  PyObject *dict = PyDict_New();
  PyDict_SetItemString(dict, "num", Py_BuildValue("i", 0));
  PyDict_SetItemString(dict, "chains", PyList_New(0));
  return dict;
}

PyObject *new_pdb() {
  PyObject *dict = PyDict_New();
  PyDict_SetItemString(dict, "name", Py_BuildValue("s", ""));
  PyDict_SetItemString(dict, "models", PyList_New(0));
  return dict;
}

PyObject *read_cif_as_pdb(PyObject *self, PyObject *args) {
  PyObject *filename_obj;
  PARSE_TUPLE(args, "O", &filename_obj);
  std::string filename = o2s(filename_obj);
  std::ifstream ifile(filename);
  jnc::bio::CifParser parser(ifile);
  std::string cat;
  std::map<std::string, std::string> values;
  cif_line_t line, old_line;

  // initialize the current atom, current residue, current chain, current model, and a pdb object
  PyObject *atom = new_atom(), *residue = new_residue(), *chain = new_chain(), *model = new_model(), *pdb = new_pdb();
  int iatom = 0;
  while (parser.next(cat, values)) {
    // if (true) {
    if (cat == "atom_site") {
      parse_cif_line(line, values);
      PyDict_SetItemString(atom, "num", Py_BuildValue("i", line.atom_num));
      PyDict_SetItemString(atom, "name", Py_BuildValue("s", line.atom_name.c_str()));
      PyList_SET_ITEM(PyDict_GetItemString(atom, "coords"), 0, Py_BuildValue("d", line.x));
      PyList_SET_ITEM(PyDict_GetItemString(atom, "coords"), 1, Py_BuildValue("d", line.y));
      PyList_SET_ITEM(PyDict_GetItemString(atom, "coords"), 2, Py_BuildValue("d", line.z));
  
      // if it's not the first atom:
      if (iatom > 0) {
        // if the new atom belongs to a different residue with the last atom:
        if (line.residue_num != old_line.residue_num || line.residue_name != old_line.residue_name || 
            line.chain_name != old_line.chain_name || line.model_num != old_line.model_num) {
          // wrap up the current residue and add to the current chain
          PyDict_SetItemString(residue, "num", Py_BuildValue("i", old_line.residue_num));
          PyDict_SetItemString(residue, "name", Py_BuildValue("s", old_line.residue_name.c_str()));
          PyList_Append(PyDict_GetItemString(chain, "residues"), residue);
          // reset the current residue
          residue = new_residue();
        }
        // if the new atom belongs to a different chain
        if (line.chain_name != old_line.chain_name) {
          // wrap up the current chain and add to the current model
          PyDict_SetItemString(chain, "name", Py_BuildValue("s", old_line.chain_name.c_str()));
          PyList_Append(PyDict_GetItemString(model, "chains"), chain);
          // reset the current chain
          chain = new_chain();
        }
        // if the new atom belongs to a different model
        if (line.model_num != old_line.model_num) {
        //   wrap up the current model and add to the pdb
          PyDict_SetItemString(model, "num", Py_BuildValue("i", old_line.model_num));
          PyList_Append(PyDict_GetItemString(pdb, "models"), model);
        //   reset the current model
          model = new_model();
        }
      }

      // add atom to the current residue and reset the current atom.
      PyList_Append(PyDict_GetItemString(residue, "atoms"), atom);
      atom = new_atom();

      old_line = std::move(line);

      iatom++;
    }
  }
  // wrap up the current residue and add to the current chain
  PyDict_SetItemString(residue, "num", Py_BuildValue("i", old_line.residue_num));
  PyDict_SetItemString(residue, "name", Py_BuildValue("s", old_line.residue_name.c_str()));
  PyList_Append(PyDict_GetItemString(chain, "residues"), residue);

  // wrap up the current chain and add to the current model
  PyDict_SetItemString(chain, "name", Py_BuildValue("s", old_line.chain_name.c_str()));
  PyList_Append(PyDict_GetItemString(model, "chains"), chain);

  // wrap up the current model and add to the pdb
  PyDict_SetItemString(model, "num", Py_BuildValue("i", old_line.model_num));
  PyList_Append(PyDict_GetItemString(pdb, "models"), model);

  return pdb;
}

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
