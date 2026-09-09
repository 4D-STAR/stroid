#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "bindings.h"

#include "stroid/config/config.h"

namespace py = pybind11;

void register_config_bindings(pybind11::module_& m) {
    py::class_<stroid::config::OptimizationMethods>(m, "OptimizationMethods")
        .def(py::init([](bool tmop, bool smoothstep) {
            return stroid::config::OptimizationMethods{tmop, smoothstep};
        }), py::arg("tmop") = false, py::arg("smoothstep") = true)
    .def_property("tmop",
        [](const stroid::config::OptimizationMethods& self) {
            return self.tmop;
        },
        [](stroid::config::OptimizationMethods& self, bool value) {
            self.tmop = value;
        }
    )
    .def_property("smoothstep",
        [](const stroid::config::OptimizationMethods& self) {
            return self.smoothstep;
        },
        [](stroid::config::OptimizationMethods& self, bool value) {
            self.smoothstep = value;
        }
    );

    py::class_<stroid::config::MeshConfig>(m, "MeshConfig")
        .def(py::init([](py::kwargs kwargs) {
            int ref_level = 4, order = 3;
            size_t continuity_order = 2, surface_bdr_id = 1, inf_bdr_id = 2, core_id = 1, envelope_id = 2, vacuum_id=3;
            bool include_external_domain = true;
            double r_core = 0.25, r_star = 1.0, flattening = 0.0, r_inf = 6.0, r_instability = 1e-14, core_steepness = 1.0;
            stroid::config::OptimizationMethods opt_method{.tmop = false, .smoothstep = true};

            return stroid::config::MeshConfig{
                .refinement_levels = kwargs.contains("refinement_levels") ? kwargs["refinement_levels"].cast<int>() : ref_level,
                .order = kwargs.contains("order") ? kwargs["order"].cast<int>() : order,
                .include_external_domain = kwargs.contains("include_external_domain") ? kwargs["include_external_domain"].cast<bool>() : include_external_domain,
                .r_core = kwargs.contains("r_core") ? kwargs["r_core"].cast<double>() : r_core,
                .r_star = kwargs.contains("r_star") ? kwargs["r_star"].cast<double>() : r_star,
                .flattening = kwargs.contains("flattening") ? kwargs["flattening"].cast<double>() : flattening,
                .r_infinity = kwargs.contains("r_infinity") ? kwargs["r_infinity"].cast<double>() : r_inf,
                .r_instability = kwargs.contains("r_instability") ? kwargs["r_instability"].cast<double>() : r_instability,
                .core_steepness = kwargs.contains("core_steepness") ? kwargs["core_steepness"].cast<double>() : core_steepness,
                .continuity_order = kwargs.contains("continuity_order") ? kwargs["continuity_order"].cast<size_t>() : continuity_order,
                .surface_bdr_id = kwargs.contains("surface_bdr_id") ? kwargs["surface_bdr_id"].cast<size_t>() : surface_bdr_id,
                .inf_bdr_id = kwargs.contains("inf_bdr_id") ? kwargs["inf_bdr_id"].cast<size_t>() : inf_bdr_id,
                .core_id = kwargs.contains("core_id") ? kwargs["core_id"].cast<size_t>() : core_id,
                .envelope_id = kwargs.contains("envelope_id") ? kwargs["envelope_id"].cast<size_t>() : envelope_id,
                .vacuum_id = kwargs.contains("vacuum_id") ? kwargs["vacuum_id"].cast<size_t>() : vacuum_id,
                .optimization_methods = kwargs.contains("optimization_methods") ? kwargs["optimization_methods"].cast<stroid::config::OptimizationMethods>() : opt_method,
                .core_mapping = kwargs.contains("core_mapping") ? kwargs["core_mapping"].cast<std::string>() : "spherified"
            };
        }))
        .def_property(
            "refinement_levels",
            [](const stroid::config::MeshConfig& self) {
                return self.refinement_levels;
            },
            [](stroid::config::MeshConfig& self, int value) {
                self.refinement_levels = value;
            }
        )
        .def_property(
            "order",
            [](const stroid::config::MeshConfig& self) {
                return self.order;
            },
            [](stroid::config::MeshConfig& self, int value) {
                self.order = value;
            }
        )
        .def_property(
            "include_external_domain",
            [](const stroid::config::MeshConfig& self) {
                return self.include_external_domain;
            },
            [](stroid::config::MeshConfig& self, bool value) {
                self.include_external_domain = value;
            }
        )
        .def_property(
            "r_core",
            [](const stroid::config::MeshConfig& self) {
                return self.r_core;
            },
            [](stroid::config::MeshConfig& self, int value) {
                self.order = value;
            }
        )
        .def_property(
            "r_star",
            [](const stroid::config::MeshConfig& self) {
                return self.r_star;
            },
            [](stroid::config::MeshConfig& self, double value) {
                self.r_star = value;
            }
        )
        .def_property(
            "flattening",
            [](const stroid::config::MeshConfig& self) {
                return self.flattening;
            },
            [](stroid::config::MeshConfig& self, double value) {
                self.flattening = value;
            }
        )
        .def_property(
            "r_infinity",
            [](const stroid::config::MeshConfig& self) {
                return self.r_infinity;
            },
            [](stroid::config::MeshConfig& self, double value) {
                self.r_infinity = value;
            }
        )
        .def_property(
            "r_instability",
            [](const stroid::config::MeshConfig& self) {
                return self.r_instability;
            },
            [](stroid::config::MeshConfig& self, double value) {
                self.r_instability = value;
            }
        )
        .def_property(
            "core_steepness",
            [](const stroid::config::MeshConfig& self) {
                return self.core_steepness;
            },
            [](stroid::config::MeshConfig& self, double value) {
                self.core_steepness = value;
            }
        )
        .def_property(
            "continuity_order",
            [](const stroid::config::MeshConfig& self) {
                return self.continuity_order;
            },
            [](stroid::config::MeshConfig& self, size_t value) {
                self.continuity_order = value;
            }
        )
        .def_property(
            "surface_bdr_id",
            [](const stroid::config::MeshConfig& self) {
                return self.surface_bdr_id;
            },
            [](stroid::config::MeshConfig& self, size_t value) {
                self.surface_bdr_id = value;
            }
        )
        .def_property(
            "inf_bdr_id",
            [](const stroid::config::MeshConfig& self) {
                return self.inf_bdr_id;
            },
            [](stroid::config::MeshConfig& self, size_t value) {
                self.inf_bdr_id = value;
            }
        )
        .def_property(
            "core_id",
            [](const stroid::config::MeshConfig& self) {
                return self.core_id;
            },
            [](stroid::config::MeshConfig& self, size_t value) {
                self.core_id = value;
            }
        )
        .def_property(
            "envelope_id",
            [](const stroid::config::MeshConfig& self) {
                return self.envelope_id;
            },
            [](stroid::config::MeshConfig& self, size_t value) {
                self.envelope_id = value;
            }
        )
        .def_property(
            "vacuum_id",
            [](const stroid::config::MeshConfig& self) {
                return self.vacuum_id;
            },
            [](stroid::config::MeshConfig& self, size_t value) {
                self.vacuum_id = value;
            }
        )
        .def_property(
            "optimization_methods",
            [](const stroid::config::MeshConfig& self) {
                return self.optimization_methods;
            },
            [](stroid::config::MeshConfig& self, stroid::config::OptimizationMethods value) {
                self.optimization_methods = value;
            }
        )
        .def_property(
            "core_mapping",
            [](const stroid::config::MeshConfig& self) {
                return self.core_mapping;
            },
            [](stroid::config::MeshConfig& self, const std::string& value) {
                if (value != "spherified" && value != "multi_block") {
                    throw std::invalid_argument("Invalid core_mapping value. Must be 'spherified' or 'multi_block'.");
                }
                self.core_mapping = value;
            }
        )
        .def("__repr__", [](const stroid::config::MeshConfig& self) {
            return stroid::config::to_string(self);
        });
}
