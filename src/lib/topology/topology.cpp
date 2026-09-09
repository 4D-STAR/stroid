#include "mfem.hpp"
#include <vector>
#include <memory>
#include <cmath>
#include <stdexcept>

#include "stroid/config/config.h"
#include "fourdst/config/config.h"

namespace stroid::topology {

    std::unique_ptr<mfem::Mesh> BuildSkeleton(const fourdst::config::Config<config::MeshConfig> & config) {
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
        mesh.FinalizeTopology();
        mesh.Finalize();
        mesh.CheckElementOrientation(true);
        mesh.CheckBdrElementOrientation(true);
        for (int i = 0; i < config->refinement_levels; ++i) {
            mesh.UniformRefinement();
        }

        if (!mesh.Conforming()) {
            std::cerr << "WARNING: Mesh has been detected to be non conforming!" << std::endl;
        }


    }

}
