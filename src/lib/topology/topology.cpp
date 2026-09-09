#include "mfem.hpp"
#include <vector>
#include <memory>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <limits>

#include "stroid/config/config.h"
#include "fourdst/config/config.h"

namespace {
    void ValidateRefinement(const stroid::config::MeshConfig& config) {
        if (config.refinement_levels.value_or(4) < 0 ||
            config.vacuum_refinement_levels.value_or(0) < 0 ||
            config.vacuum_outer_refinement_levels.value_or(0) < 0) {
            throw std::invalid_argument("Refinement levels must be non-negative.");
        }
        if (config.order.value_or(3) < 1) {
            throw std::invalid_argument("Geometry order must be at least one.");
        }

        if (!config.vacuum_refinement_levels && !config.vacuum_outer_refinement_levels) return;
        if (!config.include_external_domain.value_or(true)) {
            throw std::invalid_argument("Vacuum refinement overrides require an external domain.");
        }

        const double r_core = config.r_core.value_or(0.25);
        const double r_star = config.r_star.value_or(1.0);
        const double r_infinity = config.r_infinity.value_or(6.0);
        if (!std::isfinite(r_core) || !std::isfinite(r_star) || !std::isfinite(r_infinity) ||
            r_core <= 0.0 || r_star <= r_core || r_infinity <= r_star) {
            throw std::invalid_argument("Vacuum refinement requires finite radii with 0 < r_core < r_star < r_infinity.");
        }
        if (!std::isfinite(config.flattening.value_or(0.0)) || config.flattening.value_or(0.0) >= 1.0) {
            throw std::invalid_argument("Vacuum refinement requires finite flattening < 1.");
        }

        const auto core = config.core_id.value_or(1);
        const auto envelope = config.envelope_id.value_or(2);
        const auto vacuum = config.vacuum_id.value_or(3);
        const auto surface = config.surface_bdr_id.value_or(1);
        const auto outer = config.inf_bdr_id.value_or(2);
        const auto valid_id = [](size_t id) {
            return id > 0 && id <= static_cast<size_t>(std::numeric_limits<int>::max());
        };
        if (!valid_id(core) || !valid_id(envelope) || !valid_id(vacuum) ||
            !valid_id(surface) || !valid_id(outer) || core == envelope ||
            core == vacuum || envelope == vacuum || surface == outer) {
            throw std::invalid_argument("Vacuum refinement requires distinct positive material IDs and distinct positive boundary IDs representable as int.");
        }
    }

    void RefineReference(mfem::Mesh& mesh, const stroid::config::MeshConfig& config) {
        const int stellar_level = config.refinement_levels.value_or(4);
        const int bulk_level = config.vacuum_refinement_levels.value_or(stellar_level);
        const int outer_level = config.vacuum_outer_refinement_levels.value_or(stellar_level);
        if (!config.include_external_domain.value_or(true) ||
            (bulk_level == stellar_level && outer_level == stellar_level)) {
            for (int level = 0; level < stellar_level; ++level) mesh.UniformRefinement();
            return;
        }

        mesh.EnsureNCMesh();
        const int vacuum_id = static_cast<int>(config.vacuum_id.value_or(3));
        const double r_star = config.r_star.value_or(1.0);
        const double r_infinity = config.r_infinity.value_or(6.0);
        const double tolerance = 128.0 * std::numeric_limits<double>::epsilon() * r_infinity;

        while (true) {
            std::vector<bool> marked(static_cast<size_t>(mesh.GetNE()), false);
            for (int element = 0; element < mesh.GetNE(); ++element) {
                int target = stellar_level;
                if (mesh.GetAttribute(element) == vacuum_id) {
                    target = bulk_level;
                    double minimum = std::numeric_limits<double>::infinity();
                    double maximum = 0.0;
                    const mfem::Element* hex = mesh.GetElement(element);
                    for (int vertex = 0; vertex < hex->GetNVertices(); ++vertex) {
                        const double* position = mesh.GetVertex(hex->GetVertices()[vertex]);
                        const double radius = std::max({std::abs(position[0]), std::abs(position[1]), std::abs(position[2])});
                        minimum = std::min(minimum, radius);
                        maximum = std::max(maximum, radius);
                    }
                    if (std::abs(minimum - r_star) <= tolerance) target = std::max(target, stellar_level);
                    if (std::abs(maximum - r_infinity) <= tolerance) target = std::max(target, outer_level);
                }
                marked[static_cast<size_t>(element)] = mesh.ncmesh->GetElementDepth(element) < target;
            }

            for (int face = 0; face < mesh.GetNumFaces(); ++face) {
                const auto info = mesh.GetFaceInformation(face);
                if (!info.IsNonconformingFine() || !info.IsLocal()) continue;
                const int first = info.element[0].index;
                const int second = info.element[1].index;
                if ((mesh.GetAttribute(first) == vacuum_id) == (mesh.GetAttribute(second) == vacuum_id)) continue;
                const int coarse = mesh.ncmesh->GetElementDepth(first) < mesh.ncmesh->GetElementDepth(second) ? first : second;
                marked[static_cast<size_t>(coarse)] = true;
            }

            mfem::Array<int> refinements;
            for (int element = 0; element < mesh.GetNE(); ++element) {
                if (marked[static_cast<size_t>(element)]) refinements.Append(element);
            }
            if (refinements.Size() == 0) break;
            mesh.GeneralRefinement(refinements, 1, 1);
        }
        mesh.CheckBdrElementOrientation(true);
    }
}

namespace stroid::topology {

    std::unique_ptr<mfem::Mesh> BuildSkeleton(const fourdst::config::Config<config::MeshConfig> & config) {
        ValidateRefinement(*config);
        const std::string core_mapping = config->core_mapping.value_or("spherified");
        if (core_mapping != "spherified" && core_mapping != "multi_block") {
            throw std::invalid_argument("Unknown core_mapping: " + core_mapping);
        }

        const bool multi_block = core_mapping == "multi_block";
        const bool include_external_domain = config->include_external_domain.value_or(true);
        if (multi_block) {
            const double r_core = config->r_core.value();
            const double r_star = config->r_star.value();
            const double r_infinity = config->r_infinity.value_or(6.0);
            const double flattening = config->flattening.value();
            if (!std::isfinite(r_core) || !std::isfinite(r_star) || r_core <= 0.0 || r_star <= r_core ||
                (include_external_domain && (!std::isfinite(r_infinity) || r_infinity <= r_star))) {
                throw std::invalid_argument("multi_block requires 0 < r_core < r_star < r_infinity (when external).");
            }
            if (!std::isfinite(flattening) || flattening >= 1.0) {
                throw std::invalid_argument("multi_block requires finite flattening < 1.");
            }
        }

        const int offset = multi_block ? 8 : 0;
        int nVert = (include_external_domain ? 24 : 16) + offset;
        int nElem = (include_external_domain ? 13 : 7) + (multi_block ? 6 : 0);
        int nBev  = include_external_domain ? 12 : 6;

        auto mesh = std::make_unique<mfem::Mesh>(3, nVert, nElem, nBev, 3);

        auto add_box = [&](double scale) {
            for (const double z : {-scale, scale})
                for (const double y : {-scale, scale})
                    for (const double x : {-scale, scale})
                        mesh->AddVertex(x, y, z);
        };

        if (multi_block) {
            add_box(config->r_core.value() / 2.0);
        }
        add_box(config->r_core.value());
        add_box(config->r_star.value());
        if (include_external_domain) {
            add_box(config->r_infinity.value());
        }

        const int core_v[8] = {0, 1, 3, 2, 4, 5, 7, 6};
        mesh->AddHex(core_v, config->core_id.value());

        std::vector<std::array<int, 8>> stellar_shells = {
            {8, 9, 11, 10, 0, 1, 3, 2},
            {4, 5, 7, 6, 12, 13, 15, 14}, // +Z face
            {0, 1, 5, 4, 8, 9, 13, 12},   // -Y face
            {10, 11, 15, 14, 2, 3, 7, 6},
            {1, 3, 7, 5, 9, 11, 15, 13},  // +X face
            {0, 4, 6, 2, 8, 12, 14, 10}   // -X face
        };
        if (multi_block) {
            for (const auto & shell : stellar_shells) {
                mesh->AddHex(shell.data(), config->core_id.value());
            }
        }
        for (const auto & shell : stellar_shells) {
            auto vertices = shell;
            for (auto & vertex : vertices) vertex += offset;
            mesh->AddHex(vertices.data(), config->envelope_id.value());
        }

        if (include_external_domain) {
            std::vector<std::array<int, 8>> vacuum_shells;
            vacuum_shells.push_back({8, 9, 13, 12, 16, 17, 21, 20});
            vacuum_shells.push_back({9, 11, 15, 13, 17, 19, 23, 21});
            vacuum_shells.push_back({11, 10, 14, 15, 19, 18, 22, 23});
            vacuum_shells.push_back({10, 8, 12, 14, 18, 16, 20, 22});
            vacuum_shells.push_back({12, 13, 15, 14, 20, 21, 23, 22});
            vacuum_shells.push_back({10, 11, 9, 8, 18, 19, 17, 16});
            for (const auto & shell : vacuum_shells) {
                auto vertices = shell;
                for (auto & vertex : vertices) vertex += offset;
                mesh->AddHex(vertices.data(), config->vacuum_id.value());
            }
        }


        const int surface_bdr_quads[6][4] = {
            {12, 13, 15, 14},
            {13, 9, 11, 15},
            {9, 8, 10, 11},
            {8, 12, 14, 10},
            {8, 9, 13, 12},
            {14, 15, 11, 10}
        };

        for (const auto& bdr: surface_bdr_quads) {
            int vertices[4];
            for (int i = 0; i < 4; ++i) vertices[i] = bdr[i] + offset;
            mesh->AddBdrQuad(vertices, config->surface_bdr_id.value());
        }

        if (include_external_domain) {
            const int inf_bdr_quads[6][4] = {
                {16, 17, 21, 20},
                {17, 19, 23, 21},
                {19, 18, 22, 23},
                {18, 16, 20, 22},
                {18, 19, 17, 16},
                {20, 21, 23, 22}
            };

            for (const auto& bdr: inf_bdr_quads) {
                int vertices[4];
                for (int i = 0; i < 4; ++i) vertices[i] = bdr[i] + offset;
                mesh->AddBdrQuad(vertices, config->inf_bdr_id.value());
            }
        }

        return mesh;
    }

    // ReSharper disable once CppUseInternalLinkage
    void Finalize(mfem::Mesh& mesh, const fourdst::config::Config<config::MeshConfig> &config) {
        ValidateRefinement(*config);
        mesh.FinalizeTopology();
        mesh.Finalize();
        mesh.CheckElementOrientation(true);
        mesh.CheckBdrElementOrientation(true);
        RefineReference(mesh, *config);
    }

}
