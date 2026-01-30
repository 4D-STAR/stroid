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
     * @param config Mesh configuration (uses radii, flattening, instability radius, and core steepness).
     * @param attribute_id Element attribute ID (currently unused).
     */
    void TransformPoint(mfem::Vector& pos, const fourdst::config::Config<config::MeshConfig> &config, int attribute_id);
}