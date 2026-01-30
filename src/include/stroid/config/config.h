#pragma once

namespace stroid::config {
    /**
     * @brief Configuration parameters for stroid mesh generation.
     *
     * These values are typically loaded via
     * `fourdst::config::Config<stroid::config::MeshConfig>` from a TOML file.
     * The README shows the expected TOML layout under the `[main]` table.
     * Unspecified keys use the defaults defined here.
     */
    struct MeshConfig {
        /**
         * @brief Number of uniform refinement passes applied after topology creation.
         * @toml [main].refinement_levels
         */
        int refinement_levels = 4;
        /**
         * @brief Polynomial order for high-order elements.
         * @toml [main].order
         */
        int order = 3;
        /**
         * @brief Whether to include an external domain extending to `r_infinity`.
         * @details Currently this flag does not affect mesh generation.
         * @toml [main].include_external_domain
         */
        bool include_external_domain = false;

        /**
         * @brief Radius of the stellar core region.
         * @toml [main].r_core
         */
        double r_core = 1.5;
        /**
         * @brief Radius of the stellar surface.
         * @toml [main].r_star
         */
        double r_star = 5.0;
        /**
         * @brief Flattening factor for spheroidal shaping (0 = spherical, >0 = oblate).
         * @toml [main].flattening
         */
        double flattening = 0;

        /**
         * @brief Outer radius of the external domain when enabled.
         * @toml [main].r_infinity
         */
        double r_infinity = 6.0;

        /**
         * @brief Radius inside which transformations are skipped to avoid singularities.
         * @toml [main].r_instability
         */
        double r_instability = 1e-14;
        /**
         * @brief Controls the smoothness/steepness of the core-to-envelope transition.
         * @toml [main].core_steepness
         */
        double core_steepness = 1.0;
    };
}
