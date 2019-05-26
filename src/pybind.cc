#include "approx_hankel_transform.hpp"
#include "fast_hankel_transform.hpp"

#include <Eigen/Dense>
#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace Eigen;
using CRefCMat = const Ref<const MatrixXd>;

PYBIND11_MODULE(fhtcxx, m) {
  py::class_<FastHankelTransform>(m, "FastHankelTransform")
      .def(py::init<int, double>())
      .def("sampling", &FastHankelTransform::sampling)
      .def("set_feval", &FastHankelTransform::set_feval)
      .def("calculate", &FastHankelTransform::calculate);
  py::class_<ApproxHankelTransform>(m, "ApproxHankelTransform")
      .def(py::init<int>())
      .def("sampling", &ApproxHankelTransform::sampling)
      .def("set_feval", &ApproxHankelTransform::set_feval)
      .def("calculate", &ApproxHankelTransform::calculate);
}