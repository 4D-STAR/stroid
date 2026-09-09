#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "bindings.h"

#include "stroid/IO/mesh.h"

namespace py = pybind11;

void register_io_bindings(pybind11::module_ &m) {
    py::enum_<stroid::IO::VISUALIZATION_MODE>(m, "VISUALIZATION_MODE")
        .value("NONE", stroid::IO::VISUALIZATION_MODE::NONE)
        .value("ELEMENT_ID", stroid::IO::VISUALIZATION_MODE::ELEMENT_ID)
        .value("BOUNDARY_ELEMENT_ID", stroid::IO::VISUALIZATION_MODE::BOUNDARY_ELEMENT_ID)
        .export_values();

    m.def(
        "SaveStroidMesh",
        &stroid::IO::SaveStroidMesh,
        py::arg("mesh"),
        py::arg("filename"),
        py::arg("comment")="",
        "Save a Stroid mesh to a file."
    );
    m.def(
        "SaveMesh",
        py::overload_cast<const stroid::StroidMesh&, const std::string&>(&stroid::IO::SaveMesh),
        py::arg("mesh"),
        py::arg("filename")
    );
    m.def(
        "SaveVTU",
        py::overload_cast<const stroid::StroidMesh&, const std::string&>(&stroid::IO::SaveVTU),
        py::arg("mesh"),
        py::arg("filename")
    );
    m.def(
        "ViewMesh",
        py::overload_cast<const stroid::StroidMesh&, const std::string&, stroid::IO::VISUALIZATION_MODE, const std::string&, int, bool>(&stroid::IO::ViewMesh),
        py::arg("mesh"),
        py::arg("title")="",
        py::arg("mode")=stroid::IO::VISUALIZATION_MODE::ELEMENT_ID,
        py::arg("host")="localhost",
        py::arg("port")=19916,
        py::arg("conforming_display")=false,
        "Display the mesh in GLVis. By default, subdivide a temporary copy to avoid "
        "rendering gaps at curved hanging interfaces, preserving the source geometry "
        "and coloring. Extra display edges do not change computational DOFs. Set "
        "conforming_display=False to inspect the original element layout."
    );

    m.def(
        "VisualizeFaceValence",
        py::overload_cast<const stroid::StroidMesh&, const std::string&, int, bool>(&stroid::IO::VisualizeFaceValence),
        py::arg("mesh"),
        py::arg("host")="localhost",
        py::arg("port")=19916,
        py::arg("conforming_display")=true,
        "Display boundary-adjacent element valence: zero for untagged elements, one "
        "for surface faces, and two for internal faces (maximum if several touch an "
        "element). Values are preserved through optional display-only subdivision."
    );

    m.def(
        "ParseStroidMesh",
        [](const std::string& buf) {
            std::stringstream ss;
            ss << buf;
            auto r = stroid::IO::ParseStroidMesh(ss);
            if (!r.has_value()) {
                throw std::runtime_error("Parsing failed: " + r.error());
            }
            return std::move(r.value());
        }
    );

    m.def(
        "LoadStroidMesh",
        [](const std::string& filename) {
            auto r = stroid::IO::LoadStroidMesh(filename);
            if (!r.has_value()) {
                throw std::runtime_error("Loading " + filename + " failed: " + r.error());
            }
            return std::move(r.value());
        }
    );
}
