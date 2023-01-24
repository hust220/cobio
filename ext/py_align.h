#pragma once

#include "./py_utils.h"

// aa1 = ['A',   'R',   'N',   'D',   'C',   'Q',   'E',   'G',   'H',   'I',
// 'L',   'K',   'M',   'F',   'P',   'S',   'T',   'W',   'Y',   'V',   'O',
// 'U',   'B',   'Z',   'X',   'J',   '*'] aa2 = ['ALA', 'ARG', 'ASN', 'ASP',
// 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO',
// 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'PYL', 'SEC', 'ASX', 'GLX', 'XAA', 'XLE',
// 'TERM']

PyObject *align(PyObject *self, PyObject *args);

PyObject *map_atoms(PyObject *self, PyObject *args);

PyObject *rmsd(PyObject *self, PyObject *args);

PyObject *suppos(PyObject *self, PyObject *args);

// Apply superposition to a list of residues
PyObject *sp_apply(PyObject *self, PyObject *args);





