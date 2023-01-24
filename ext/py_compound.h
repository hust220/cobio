#pragma once

#include "./py_utils.h"

extern PyTypeObject PocketType;

PyObject *compound_grow(PyObject *self, PyObject *args, PyObject *kwargs);

