#pragma once

#include "mfem.hpp"
#include "stroid/config/config.h"
#include "fourdst/config/config.h"

namespace stroid::topology {
    /**
     * @brief Apply an equiangular (gnomonic) projection to a point on a cube.
     * @param pos Position vector updated in-place.
     */
    void ApplyEquiangular(mfem::Vector& pos);

    /**
     * @brief Apply spheroidal flattening along the Z axis.
     * @param pos Position vector updated in-place.
     * @param config Mesh configuration (uses `flattening`).
     */
    void ApplySpheroidal(mfem::Vector& pos, const fourdst::config::Config<config::MeshConfig> &config);

    /**
     * @brief Apply Kelvin transform outside the stellar radius.
     * @param pos Position vector updated in-place.
     * @param config Mesh configuration (uses `r_star` and `r_infinity`).
     */
    void ApplyKelvin(mfem::Vector& pos, const fourdst::config::Config<config::MeshConfig> &config);

    /**
     * @brief Map a point from the initial block topology to the curvilinear domain.
     * @param pos Position vector updated in-place.
     * @param config Mesh configuration (uses radii, flattening, and `core_mapping`).
     * The `multi_block` strategy requires the matching skeleton from BuildSkeleton;
     * changing only the mapping on a legacy core element is not supported.
     * @param attribute_id Element attribute ID (currently unused).
     */
    void TransformPoint(mfem::Vector& pos, const fourdst::config::Config<config::MeshConfig> &config, int attribute_id);

    /**
     * @brief Compute the compactification coordinate for a point in the curvilinear domain. This ranges from 0-1 with 0 at the stellar surface and 1 at the compactified infinity.
     * @param logical_position Logical position of the point in the curvilinear domain.
     * @param attribute Element attribute ID (used to determine the exterior coordinate).
     * @param config Mesh configuration (uses radii and flattening).
     * @return Compactification coordinate ranging from 0 (stellar surface) to 1 (compactified infinity).
     */
    double ComputeExteriorCoordinate(
        const mfem::Vector& logical_position,
        int attribute,
        const fourdst::config::Config<config::MeshConfig>& config
    );
}
