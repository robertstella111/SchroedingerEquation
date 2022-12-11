#include "io.hpp"

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>
#include <memory>
#include <complex>
#include <array>
#include <cmath>
namespace py = pybind11;

PYBIND11_MODULE(Schroedinger, m) {

  py::class_<Schroedinger1D>(m, "Schroedinger1D")
      .def(py::init<>())
      .def("solve", &Schroedinger1D::solve, "solve")
      .def("getEigenvektor", &Schroedinger1D::getEigenvektorSort, "getEigenvektor")
      .def("getXKoord", &Schroedinger1D::getXKoord, "getXKoord")
      .def("setNumPoints", &Schroedinger1D::setNumPoints, "setLength  ")
      .def("getEigenValues", &Schroedinger1D::getEigenValuesSort, "setLength  ")
      .def("setDomain", &Schroedinger1D::setDomain, "setLength  ");





      

}
