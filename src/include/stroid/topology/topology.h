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
     * @brief Finalize topology, validate orientation, and apply uniform refinement.
     * @param mesh Mesh to finalize in-place.
     * @param config Mesh configuration (uses `refinement_levels`).
     */
    void Finalize(mfem::Mesh& mesh, const fourdst::config::Config<config::MeshConfig> &config);
}
