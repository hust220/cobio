#include <pybind11/pybind11.h>
#include "bio.h"

namespace py = pybind11;

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

class Array3d {
public:
    std::array<int, 3> shape;
    std::array<int, 3> strides;
    double *data;

    Array3d(const std::array<int, 3> &shape_) : shape(shape_) {
        data = new double[shape[0]*shape[1]*shape[2]];
        strides[0] = shape[1]*shape[2];
        strides[1] = shape[2];
        strides[2] = 1;
    }

    double &operator ()(int i, int j, int k) {
        return data[i*strides[0]+j*strides[1]+k*strides[2]];
    }
};

int add(int i, int j) {
    return i + j;
}

Array3d extract(const std::string &recfile, const std::string &ligfile) {
    auto rec = jnc::bio::read_pdb(recfile);
    auto lig = jnc::bio::read_mol2s(ligfile)[0];
    auto rs = rec[0].presidues();

    int nresidues = rs.size();
    int lig_size = lig.atoms.size();

    std::array<double, 3> c {0, 0, 0};
    for (auto && atom : lig.atoms) {
        c[0] += atom.x;
        c[1] += atom.y; 
        c[2] += atom.z;
    }
    for (int i = 0; i < 3; i++) c[i] /= double(lig_size);

    int bins = 64;
    double size = 32;
    int bin_size = bins / size;
    Array3d m({bins, bins, bins});
    std::array<double, 3> ll { c[0]-size/2, c[1]-size/2, c[2]-size/2 }; // lower limit
    for (auto && atom : lig.atoms) {
        int dx = int((atom.x-ll[0]) / bin_size);
        int dy = int((atom.y-ll[1]) / bin_size);
        int dz = int((atom.z-ll[2]) / bin_size);
        if (dx < bins && dy < bins && dz < bins) {
            m(dx, dy, dz) = 1;
        }
    }
    for (auto && pr : rs) {
        for (auto && atom : *pr) {
            int dx = int((atom[0]-ll[0]) / bin_size);
            int dy = int((atom[1]-ll[1]) / bin_size);
            int dz = int((atom[2]-ll[2]) / bin_size);
            if (dx < bins && dy < bins && dz < bins) {
                m(dx, dy, dz) = 1;
            }
        }
    }
    return m;
}

namespace py = pybind11;

PYBIND11_MODULE(mlpocket, m) {
    py::class_<Array3d>(m, "Array3d", py::buffer_protocol())
        .def_buffer([](Array3d &m) -> py::buffer_info {
            return py::buffer_info(
                m.data,                               /* Pointer to buffer */
                sizeof(double),                          /* Size of one scalar */
                py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
                3,                                      /* Number of dimensions */
                { m.shape[0], m.shape[1], m.shape[2] },                 /* Buffer dimensions */
                { sizeof(double)*m.strides[0], sizeof(double)*m.strides[1], sizeof(double)*m.strides[2] }
                );
        });

    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: mlpocket

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("extract", &extract, R"pbdoc(
        Read PDB file

        Some other explanation about the add function.
    )pbdoc");

    m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
