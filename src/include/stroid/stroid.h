#pragma once

#include "stroid/config/config.h"
#include "stroid/topology/topology.h"
#include "stroid/topology/mapping.h"
#include "stroid/topology/curvilinear.h"
#include "stroid/topology/optimize.h"
#include "stroid/utils/mesh_utils.h"
#include "stroid/IO/mesh.h"
#include "stroid/utils/types.h"
#include "stroid/refinement/uniform.h"
#include "stroid/utils/mesh_stats.h"
#include "stroid/version.h"

/**
 * @namespace stroid
 * @brief Public API surface for the stroid mesh generation library.
 *
 * This namespace aggregates configuration, topology construction, mapping, and
 * mesh utilities for building curvilinear multi-block meshes.
 *
 * @par Example: build a mesh programmatically
 * @code
 * #include <memory>
 * #include "mfem.hpp"
 * #include "fourdst/config/config.h"
 * #include "stroid/stroid.h"
 *
 * int main() {
 *     fourdst::config::Config<stroid::config::MeshConfig> cfg;
 *     auto mesh = stroid::topology::BuildSkeleton(cfg);
 *     stroid::topology::Finalize(*mesh, cfg);
 *     stroid::topology::PromoteToHighOrder(*mesh, cfg);
 *     stroid::topology::ProjectMesh(*mesh, cfg);
 *     stroid::IO::SaveVTU(*mesh, "stroid");
 * }
 * @endcode
 *
 * @par Example: tweak config then visualize
 * @code
 * fourdst::config::Config<stroid::config::MeshConfig> cfg;
 * cfg->order = 4;
 * cfg->flattening = 0.1;
 * auto mesh = stroid::topology::BuildSkeleton(cfg);
 * stroid::topology::Finalize(*mesh, cfg);
 * stroid::topology::PromoteToHighOrder(*mesh, cfg);
 * stroid::topology::ProjectMesh(*mesh, cfg);
 * stroid::IO::ViewMesh(*mesh, "Stroid Mesh", stroid::IO::VISUALIZATION_MODE::ELEMENT_ID, "localhost", 19916);
 * @endcode
 */
namespace stroid {
    inline StroidMesh GenerateMesh(const fourdst::config::Config<stroid::config::MeshConfig>& input) {
        fourdst::config::Config<config::MeshConfig> cfg;
        cfg.mutate([&input](config::MeshConfig& value) {
            value = config::ResolveDefaults(*input);
        });
        StroidMesh sm;
        sm.type = MFEM_MESH_TYPE::SERIAL;
        sm.config = *cfg;
        auto reference = stroid::topology::BuildSkeleton(cfg);
        stroid::topology::Finalize(*reference, cfg);
        sm.refinement_levels = cfg->refinement_levels.value_or(0);

        sm.reference_mesh = std::move(reference);
        sm.mesh = utils::BuildProjected(*sm.reference_mesh, cfg);
        stroid::topology::OptimizeMesh(*sm.mesh, cfg);
        sm.exterior_coordinate = stroid::topology::BuildExteriorCoordinate(*sm.mesh, *sm.reference_mesh, cfg);
        return sm;
    }
    inline StroidMesh GenerateMesh(const stroid::config::MeshConfig& config) {
        fourdst::config::Config<config::MeshConfig> cfg;
        auto Mutator = [&config](config::MeshConfig& orig) {
            orig = config;
        };

        cfg.mutate(Mutator);
        return GenerateMesh(cfg);
    }
    inline StroidMesh GenerateMesh(const std::string& filename) {
        fourdst::config::Config<stroid::config::MeshConfig> config;
        config.load(filename);
        return GenerateMesh(config);
    }
}


/**
 * @namespace stroid::config
 * @brief Configuration types and defaults for mesh generation.
 */
namespace stroid::config {}

/**
 * @namespace stroid::topology
 * @brief Topology construction and curvilinear mapping utilities.
 */
namespace stroid::topology {}

/**
 * @namespace stroid::IO
 * @brief Mesh serialization and visualization helpers.
 */
namespace stroid::IO {}

/**
 * @namespace stroid::utils
 * @brief Mesh inspection and validation utilities.
 */
namespace stroid::utils {}
