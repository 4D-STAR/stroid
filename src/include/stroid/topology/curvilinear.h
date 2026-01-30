#pragma once

#include "mfem.hpp"
#include "stroid/config/config.h"
#include "fourdst/config/config.h"

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
}