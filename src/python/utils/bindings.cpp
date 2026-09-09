#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "bindings.h"

#include "stroid/utils/types.h"
#include "stroid/utils/mesh_stats.h"
#include "stroid/utils/mesh_utils.h"

namespace py = pybind11;

void register_stats_bindings(pybind11::module_ &m) {
    auto statsMod = m.def_submodule("stats", "Stats Bindings");
    py::enum_<stroid::stats::MeshStatFeatures>(statsMod, "MeshStatFeatures", py::arithmetic())
     .value("NONE", stroid::stats::MeshStatFeatures::NONE)
     .value("RADIUS", stroid::stats::MeshStatFeatures::RADIUS)
     .value("AXES", stroid::stats::MeshStatFeatures::AXES)
     .value("ELLIPTICITY", stroid::stats::MeshStatFeatures::ELLIPTICITY)
     .value("BOWING", stroid::stats::MeshStatFeatures::BOWING)
     .value("CONFORMITY", stroid::stats::MeshStatFeatures::CONFORMITY)
     .value("JACOBIAN", stroid::stats::MeshStatFeatures::JACOBIAN)
     .value("VOLUME_AREA", stroid::stats::MeshStatFeatures::VOLUME_AREA)
     .value("ELEMENT_COUNT", stroid::stats::MeshStatFeatures::ELEMENT_COUNT)
     .value("MESH_SIZE", stroid::stats::MeshStatFeatures::MESH_SIZE)
     .value("OUTER_BOUNDS", stroid::stats::MeshStatFeatures::OUTER_BOUNDS)
     .value("CENTROID", stroid::stats::MeshStatFeatures::CENTROID)
     .value("CONFIG_META", stroid::stats::MeshStatFeatures::CONFIG_META)
     .value("BOUNDING_BOX", stroid::stats::MeshStatFeatures::BOUNDING_BOX)
     .value("REFINEMENT", stroid::stats::MeshStatFeatures::REFINEMENT)
     .export_values()
     .def("__or__", [](stroid::stats::MeshStatFeatures lhs, stroid::stats::MeshStatFeatures rhs) {
         return lhs | rhs;
     }, py::is_operator())
     .def("__and__", [](stroid::stats::MeshStatFeatures lhs, stroid::stats::MeshStatFeatures rhs) {
         return lhs & rhs;
     }, py::is_operator());

    py::class_<stroid::stats::RadiusStats>(statsMod, "RadiusStats")
    .def_readonly("min", &stroid::stats::RadiusStats::min)
    .def_readonly("max", &stroid::stats::RadiusStats::max)
    .def_readonly("mean", &stroid::stats::RadiusStats::mean)
    .def_readonly("stddev", &stroid::stats::RadiusStats::stddev)
    .def_readonly("n_samples", &stroid::stats::RadiusStats::n_samples);

    py::class_<stroid::stats::AxisStats>(statsMod, "AxisStats")
    .def_readonly("semi_major", &stroid::stats::AxisStats::semi_major)
    .def_readonly("semi_minor", &stroid::stats::AxisStats::semi_minor);

    py::class_<stroid::stats::EllipticityStats>(statsMod, "EllipticityStats")
    .def_readonly("flattening", &stroid::stats::EllipticityStats::flattening)
    .def_readonly("polar_equatorial", &stroid::stats::EllipticityStats::polar_equatorial)
    .def_readonly("radius_uniformity", &stroid::stats::EllipticityStats::radius_uniformity);

    py::class_<stroid::stats::BowingStats>(statsMod, "BowingStats")
    .def_readonly("max_inward", &stroid::stats::BowingStats::max_inward)
    .def_readonly("max_outward", &stroid::stats::BowingStats::max_outward)
    .def_readonly("rms", &stroid::stats::BowingStats::rms)
    .def_readonly("worst_at_radius", &stroid::stats::BowingStats::worst_at_radius);

    py::class_<stroid::stats::ConformityStats>(statsMod, "ConformityStats")
    .def_readonly("conforming", &stroid::stats::ConformityStats::conforming)
    .def_readonly("hierarchy_enabled", &stroid::stats::ConformityStats::hierarchy_enabled)
    .def_readonly("n_nonconforming_faces", &stroid::stats::ConformityStats::n_nonconforming_faces);

    py::class_<stroid::stats::RegionRefinementStats>(statsMod, "RegionRefinementStats")
    .def_readonly("min_depth", &stroid::stats::RegionRefinementStats::min_depth)
    .def_readonly("max_depth", &stroid::stats::RegionRefinementStats::max_depth);

    py::class_<stroid::stats::RefinementStats>(statsMod, "RefinementStats")
    .def_readonly("all", &stroid::stats::RefinementStats::all)
    .def_readonly("core", &stroid::stats::RefinementStats::core)
    .def_readonly("envelope", &stroid::stats::RefinementStats::envelope)
    .def_readonly("vacuum", &stroid::stats::RefinementStats::vacuum)
    .def_readonly("geometry_dofs", &stroid::stats::RefinementStats::geometry_dofs)
    .def_readonly("geometry_true_dofs", &stroid::stats::RefinementStats::geometry_true_dofs);

    py::class_<stroid::stats::JacobianStats>(statsMod, "JacobianStats")
    .def_readonly("detJ_min", &stroid::stats::JacobianStats::detJ_min)
    .def_readonly("detJ_max", &stroid::stats::JacobianStats::detJ_max)
    .def_readonly("n_flipped", &stroid::stats::JacobianStats::n_flipped)
    .def_readonly("min_detJ_ratio", &stroid::stats::JacobianStats::min_detJ_ratio)
    .def_readonly("worst_ratio_at_radius", &stroid::stats::JacobianStats::worst_ratio_at_radius)
    .def_readonly("detJ_min_at_radius", &stroid::stats::JacobianStats::detJ_min_at_radius)
    .def_readonly("n_elements", &stroid::stats::JacobianStats::n_elements);

    py::class_<stroid::stats::VolumeAreaStats>(statsMod, "VolumeAreaStats")
    .def_readonly("stellar_volume", &stroid::stats::VolumeAreaStats::stellar_volume)
    .def_readonly("surface_area", &stroid::stats::VolumeAreaStats::surface_area)
    .def_readonly("analytic_volume", &stroid::stats::VolumeAreaStats::analytic_volume)
    .def_readonly("analytic_area", &stroid::stats::VolumeAreaStats::analytic_area);

    py::class_<stroid::stats::ElementCounts>(statsMod, "ElementCounts")
    .def_readonly("total", &stroid::stats::ElementCounts::total)
    .def_readonly("core", &stroid::stats::ElementCounts::core)
    .def_readonly("envelope", &stroid::stats::ElementCounts::envelope)
    .def_readonly("vacuum", &stroid::stats::ElementCounts::vacuum)
    .def_readonly("other", &stroid::stats::ElementCounts::other)
    .def_readonly("n_vertices", &stroid::stats::ElementCounts::n_vertices);

    py::class_<stroid::stats::MeshSizeStats>(statsMod, "MeshSizeStats")
    .def_readonly("h_min", &stroid::stats::MeshSizeStats::h_min)
    .def_readonly("h_max", &stroid::stats::MeshSizeStats::h_max)
    .def_readonly("h_mean", &stroid::stats::MeshSizeStats::h_mean)
    .def_readonly("h_stddev", &stroid::stats::MeshSizeStats::h_stddev);

    py::class_<stroid::stats::OuterBoundsStats>(statsMod, "OuterBoundsStats")
    .def_readonly("min", &stroid::stats::OuterBoundsStats::min)
    .def_readonly("max", &stroid::stats::OuterBoundsStats::max)
    .def_readonly("mean", &stroid::stats::OuterBoundsStats::mean)
    .def_readonly("n_samples", &stroid::stats::OuterBoundsStats::n_samples);

    py::class_<stroid::stats::CentroidStats>(statsMod, "CentroidStats")
    .def_readonly("x", &stroid::stats::CentroidStats::x)
    .def_readonly("y", &stroid::stats::CentroidStats::y)
    .def_readonly("z", &stroid::stats::CentroidStats::z)
    .def_readonly("offset", &stroid::stats::CentroidStats::offset);

    py::class_<stroid::stats::ConfigMeta>(statsMod, "ConfigMeta")
    .def_readonly("r_core", &stroid::stats::ConfigMeta::r_core)
    .def_readonly("r_star", &stroid::stats::ConfigMeta::r_star)
    .def_readonly("flattening", &stroid::stats::ConfigMeta::flattening)
    .def_readonly("r_infinity", &stroid::stats::ConfigMeta::r_infinity)
    .def_readonly("geom_order", &stroid::stats::ConfigMeta::geom_order)
    .def_readonly("refinement_levels", &stroid::stats::ConfigMeta::refinement_levels)
    .def_readonly("has_external_domain", &stroid::stats::ConfigMeta::has_external_domain);

    py::class_<stroid::stats::BoundingBox>(statsMod, "BoundingBox")
    .def_readonly("xMin", &stroid::stats::BoundingBox::xMin)
    .def_readonly("xMax", &stroid::stats::BoundingBox::xMax)
    .def_readonly("yMin", &stroid::stats::BoundingBox::yMin)
    .def_readonly("yMax", &stroid::stats::BoundingBox::yMax)
    .def_readonly("zMin", &stroid::stats::BoundingBox::zMin)
    .def_readonly("zMax", &stroid::stats::BoundingBox::zMax)
    .def_readonly("valid", &stroid::stats::BoundingBox::valid)
    .def("dx", &stroid::stats::BoundingBox::dx)
    .def("dy", &stroid::stats::BoundingBox::dy)
    .def("dz", &stroid::stats::BoundingBox::dz)
    .def("diag", &stroid::stats::BoundingBox::diag);

    py::class_<stroid::stats::BoundingBoxStats>(statsMod, "BoundingBoxStats")
    .def_readonly("core", &stroid::stats::BoundingBoxStats::core)
    .def_readonly("star", &stroid::stats::BoundingBoxStats::star)
    .def_readonly("vacuum", &stroid::stats::BoundingBoxStats::vacuum);

    py::class_<stroid::stats::MeshStats>(statsMod, "MeshStats")
    .def_readonly("computed", &stroid::stats::MeshStats::computed)
    .def_readonly("radius", &stroid::stats::MeshStats::radius)
    .def_readonly("axes", &stroid::stats::MeshStats::axes)
    .def_readonly("ellipticity", &stroid::stats::MeshStats::ellipticity)
    .def_readonly("bowing", &stroid::stats::MeshStats::bowing)
    .def_readonly("conformity", &stroid::stats::MeshStats::conformity)
    .def_readonly("refinement", &stroid::stats::MeshStats::refinement)
    .def_readonly("jacobian", &stroid::stats::MeshStats::jacobian)
    .def_readonly("jacobian_stellar", &stroid::stats::MeshStats::jacobian_stellar)
    .def_readonly("jacobian_vacuum", &stroid::stats::MeshStats::jacobian_vacuum)
    .def_readonly("volume", &stroid::stats::MeshStats::volume)
    .def_readonly("element_counts", &stroid::stats::MeshStats::element_counts)
    .def_readonly("mesh_size", &stroid::stats::MeshStats::mesh_size)
    .def_readonly("outer_bounds", &stroid::stats::MeshStats::outer_bounds)
    .def_readonly("centroid", &stroid::stats::MeshStats::centroid)
    .def_readonly("config_meta", &stroid::stats::MeshStats::config_meta)
    .def_readonly("bounding_box", &stroid::stats::MeshStats::bounding_box)
    .def_readonly("warnings", &stroid::stats::MeshStats::warnings)
    .def_readonly("errors", &stroid::stats::MeshStats::errors)
    .def("__repr__", [](const stroid::stats::MeshStats& self) {
       return stroid::stats::to_string(self);
    });

    statsMod.attr("MESH_STAT_DEFAULT") = stroid::stats::MESH_STAT_DEFAULT;
    statsMod.attr("MESH_STAT_ALL") = stroid::stats::MESH_STAT_ALL;

    statsMod.def(
        "ComputeMeshStats",
        &stroid::stats::ComputeMeshStats,
        py::arg("mesh"),
        py::arg("features") = stroid::stats::MESH_STAT_DEFAULT,
        py::arg("sample_order")=-1
    );
}

void register_type_bindings(py::module_ &m) {
    py::enum_<stroid::MFEM_MESH_TYPE>(m, "MFEM_MESH_TYPE")
    .value("SERIAL", stroid::MFEM_MESH_TYPE::SERIAL)
    .value("PARALLEL", stroid::MFEM_MESH_TYPE::PARALLEL)
    .export_values();

    py::class_<stroid::StroidMesh>(m, "StroidMesh")
    .def_property_readonly("type", [](const stroid::StroidMesh& self) {
        return (self.type == stroid::MFEM_MESH_TYPE::SERIAL) ? "SERIAL" : "PARALLEL";
    })
    .def_readonly("config", &stroid::StroidMesh::config)
    .def_readonly("refinement_levels", &stroid::StroidMesh::refinement_levels)
    .def("has_mesh", [](const stroid::StroidMesh& self) {
        return self.mesh != nullptr;
    })
    .def("has_rmesh", [](const stroid::StroidMesh& self) {
        return self.reference_mesh != nullptr;
    })
    .def("mesh_stats", [](const stroid::StroidMesh& self, bool use_ref_mesh) {
        auto result = self.mesh_stats(use_ref_mesh);
        if (!result.has_value()) {
            throw std::runtime_error(result.error());
        }
        return result.value();
    }, py::arg("use_ref_mesh") = false)
    .def("__repr__", [](const stroid::StroidMesh& self) {
       return std::format("<StroidMesh [{}]: NE: {}, NV: {}>", (self.type == stroid::MFEM_MESH_TYPE::SERIAL) ? "SERIAL" : "PARALLEL", self.mesh->GetNE(), self.mesh->GetNV());
    });
}

void register_utils_bindings(pybind11::module_ &m) {
    register_type_bindings(m);
    register_stats_bindings(m);
}
