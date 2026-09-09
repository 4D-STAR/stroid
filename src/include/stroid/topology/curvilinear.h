#pragma once

#include "mfem.hpp"
#include "stroid/config/config.h"
#include "fourdst/config/config.h"
#include "stroid/utils/types.h"

namespace stroid::topology {
    /**
     * @brief Promote a mesh to high-order by attaching an H1 nodal finite element space.
     * @param mesh Mesh to update in-place.
     * @param config Mesh configuration (uses `order`).
     */
    void PromoteToHighOrder(mfem::Mesh& mesh, const fourdst::config::Config<config::MeshConfig> &config);
    /**
     * @brief Project high-order mesh nodes using the configured curvilinear mapping.
     * @details Requires nodes to be present (call PromoteToHighOrder first).
     * @param mesh Mesh to update in-place.
     * @param config Mesh configuration (uses radii, flattening, and mapping parameters).
     */
    void ProjectMesh(mfem::Mesh& mesh, const fourdst::config::Config<config::MeshConfig> &config);

    /**
     * @brief Build a scalar grid function representing the compactification coordinate for a mesh. This ranges from 0-1 with 0 at the stellar surface and 1 at the compactified infinity.
     * @param mesh Reference to the underlying serial MFEM mesh which has been promoted to high-order and projected into the curvilinear domain.
     * @param reference_mesh reference to the underlying serial which has not been promoted to high-order or projected into the curvilinear domain. This is used to compute the compactification coordinate.
     * @param config Config file
     * @return Unique pointer to a scalar mesh field representing the compactification coordinate.
     */
    std::unique_ptr<ScalarMeshField> BuildExteriorCoordinate(
        mfem::Mesh& mesh,
        mfem::Mesh& reference_mesh,
        const fourdst::config::Config<config::MeshConfig>& config
    );
}