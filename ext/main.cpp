#include "./py_align.h"
#include "./py_compound.h"
#include "./py_pdb.h"
#include "./py_utils.h"

static PyMethodDef cobio_methods[] = {
    // {"add", add, METH_VARARGS, "Add"},
    //    {"newatom", newatom, METH_VARARGS, "New Atom"},
    {"suppos", suppos, METH_VARARGS, "Superposion"},
    {"sp_apply", sp_apply, METH_VARARGS, "Apply Superposion"},
    {"align", align, METH_VARARGS, "Apply Superposion"},
    {"rmsd", rmsd, METH_VARARGS, "RMSD"},
    {"map_atoms", map_atoms, METH_VARARGS, "Map Atoms"},
    {"read_cif_as_pdb", read_cif_as_pdb, METH_VARARGS, "Read cif file as pdb dictionary"},
    {"compound_grow", (PyCFunction)compound_grow, METH_VARARGS | METH_KEYWORDS, "Grow a compound atom by atom"},
    //    {"Amorphize", amorphize, METH_VARARGS, "Amorphize"},
    //    {"aa321", aa321, METH_VARARGS, "aa321"},
    //    {"aa123", aa123, METH_VARARGS, "aa123"},
    {NULL, NULL, 0, NULL} // sentinel
};

static PyModuleDef cobio_module = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "_cobio",
    .m_doc = "Jian's Python Library",
    .m_size = -1,
    .m_methods = cobio_methods,
};

PyMODINIT_FUNC PyInit__cobio(void) {
  // return PyModule_Create(&cobio_module);

  PyObject *m;

  if (PyType_Ready(&PocketType) < 0)
    return NULL;

  if (PyType_Ready(&PdbType) < 0)
    return NULL;

  m = PyModule_Create(&cobio_module);
  if (m == NULL)
    return NULL;

  Py_INCREF(&PocketType);
  if (PyModule_AddObject(m, "Pocket", (PyObject *)&PocketType) < 0) {
    Py_DECREF(&PocketType);
    Py_DECREF(m);
    return NULL;
  }

  Py_INCREF(&PdbType);
  if (PyModule_AddObject(m, "Pdb", (PyObject *)&PdbType) < 0) {
    Py_DECREF(&PdbType);
    Py_DECREF(m);
    return NULL;
  }

  return m;
}
