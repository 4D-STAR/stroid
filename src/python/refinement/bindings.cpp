#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "bindings.h"

#include "stroid/refinement/uniform.h"

namespace py = pybind11;

void register_refinement_bindings(pybind11::module_ &m) {
    m.def("UniformRefinement", &stroid::refinement::UniformRefinement, py::arg("mesh"), py::arg("levels"), "Refine every current leaf by the requested additional levels, preserving hanging-node constraints and rebuilding high-order geometry. Initial configuration targets remain unchanged.");
}
