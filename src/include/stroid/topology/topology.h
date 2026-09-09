#pragma once
#include "mfem.hpp"
#include "stroid/config/config.h"
#include "fourdst/config/config.h"

#include <memory>

namespace stroid::topology {
    /**
     * @brief Build the initial multi-block mesh topology for the star model.
     * @param config Mesh configuration (uses radii, domain flags, and `core_mapping`).
     * The legacy `spherified` core uses one block; `multi_block` uses an
     * inner Cartesian block and six core transition blocks.
     * @return Newly allocated mesh skeleton (not yet refined or curved).
     */
    std::unique_ptr<mfem::Mesh> BuildSkeleton(const fourdst::config::Config<config::MeshConfig> & config);
    /**
     * @brief Finalize topology, validate orientation, and apply the configured refinement policy.
     * @param mesh Mesh to finalize in-place.
     * @param config Mesh configuration. Vacuum refinement overrides enable a
     * balanced hierarchy with fine layers at both vacuum boundaries and a
     * conforming stellar interface. Without overrides, refinement is uniform.
     */
    void Finalize(mfem::Mesh& mesh, const fourdst::config::Config<config::MeshConfig> &config);
}
