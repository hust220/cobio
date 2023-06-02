#pragma once

#include "./py_utils.h"

// aa1 = ['A',   'R',   'N',   'D',   'C',   'Q',   'E',   'G',   'H',   'I',
// 'L',   'K',   'M',   'F',   'P',   'S',   'T',   'W',   'Y',   'V',   'O',
// 'U',   'B',   'Z',   'X',   'J',   '*'] aa2 = ['ALA', 'ARG', 'ASN', 'ASP',
// 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO',
// 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'PYL', 'SEC', 'ASX', 'GLX', 'XAA', 'XLE',
// 'TERM']

extern PyTypeObject PdbType;


PyObject *read_cif_as_pdb(PyObject *self, PyObject *args);


