#pragma once

#define PY_SSIZE_T_CLEAN

#include "Python.h"
#include "./jnc_core.h"
#include "./align.h"
#include <array>
#include <numeric>
#include <stdio.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

template <typename T> PyObject *to_obj(T);
template <> inline PyObject *to_obj(int n) { return PyLong_FromLong(n); }
template <> inline PyObject *to_obj(double n) { return PyFloat_FromDouble(n); }
template <> inline PyObject *to_obj(const std::string &str) { return PyUnicode_FromString(str.c_str()); }

template <typename T> T obj_to(PyObject *);
template <> inline int obj_to<int>(PyObject *o) { return int(PyLong_AS_LONG(o)); }
template <> inline double obj_to<double>(PyObject *o) { return PyFloat_AS_DOUBLE(o); }
template <> inline std::string obj_to<std::string>(PyObject *o) { return PyBytes_AS_STRING(o); }
template <> inline PyObject *obj_to<PyObject *>(PyObject *o) { return o; }

template <typename Output_, typename Input_> Output_ obj_get(PyObject *obj, const Input_ &key) {
  auto key_object = to_obj(key);
  auto item_object = PyObject_GetItem(obj, key_object);
  Py_DECREF(key_object);

  auto item = obj_to<Output_>(item_object);
  Py_DECREF(item_object);
  return item;
}

struct MatViewer {
  PyObject *obj;
  MatViewer(PyObject *o) : obj(o) {}
  double operator()(int i, int j) const { return PyFloat_AS_DOUBLE(PyList_GET_ITEM(PyList_GET_ITEM(obj, i), j)); }
};

inline double at(PyObject *vec, int i) { return PyFloat_AS_DOUBLE(PyList_GET_ITEM(vec, i)); }

inline void assign(PyObject *vec, int i, double v) { PyList_SET_ITEM(vec, i, PyFloat_FromDouble(v)); }

inline char *o2s(PyObject *o) {
  if (PyUnicode_Check(o)) {
    o = PyUnicode_AsUTF8String(o);
  }
  if (!PyBytes_Check(o)) {
    PyErr_SetString(PyExc_RuntimeError, "Not a string!");
  }
  return PyBytes_AS_STRING(o);
}

inline int o2i(PyObject *o) { return int(PyLong_AsLong(PyNumber_Long(o))); }

inline double o2d(PyObject *o) { return PyFloat_AS_DOUBLE(PyNumber_Float(o)); }

template <typename Array_> PyObject *a2o(const Array_ &c, int size) {
  PyObject *obj = PyList_New(size);
  for (int i = 0; i < size; i++) {
    PyList_SET_ITEM(obj, i, PyFloat_FromDouble(c[i]));
  }
  return obj;
}

inline PyObject *attr(PyObject *o, const char *key) { return PyObject_GetAttrString(o, key); }

#define Err(...) PyErr_SetString(PyExc_TypeError, __VA_ARGS__)

#define PARSE_TUPLE(...)                                                                                               \
  if (!PyArg_ParseTuple(__VA_ARGS__)) {                                                                                \
    Err("parameter type error");                                                                                       \
    return NULL;                                                                                                       \
  }

inline PyObject *str2obj(const char *str) {
  //    return PyUnicode_DecodeUTF8(PyBytes_FromString(str));
  return Py_BuildValue("s", str);
}

template <typename _Mat> void rotate(PyObject *o, const _Mat &mat) {
  double _x = at(o, 0);
  double _y = at(o, 1);
  double _z = at(o, 2);
  assign(o, 0, _x * mat(0, 0) + _y * mat(1, 0) + _z * mat(2, 0));
  assign(o, 1, _x * mat(0, 1) + _y * mat(1, 1) + _z * mat(2, 1));
  assign(o, 2, _x * mat(0, 2) + _y * mat(1, 2) + _z * mat(2, 2));
}