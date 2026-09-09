#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <format>
#include <sstream>

namespace stroid::config {

    struct OptimizationMethods {
        std::optional<bool> tmop{false};
        std::optional<bool> smoothstep{true};
    };

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
         * @section toml
         * - [main].refinement_levels
         */
        std::optional<int> refinement_levels = 4;
        /**
         * @brief Polynomial order for high-order elements.
         * @section toml
         * - [main].order
         */
        std::optional<int> order = 3;
        /**
         * @brief Whether to include an external domain extending to `r_infinity`.
         * @section toml
         * - [main].include_external_domain
         */
        std::optional<bool> include_external_domain = true;

        /**
         * @brief Radius of the stellar core region.
         * @section toml
         * - [main].r_core
         */
        std::optional<double> r_core = 0.25;
        /**
         * @brief Radius of the stellar surface.
         * @section toml
         * - [main].r_star
         */
        std::optional<double> r_star = 1.0;
        /**
         * @brief Flattening factor for spheroidal shaping (0 = spherical, >0 = oblate).
         * @section toml
         * - [main].flattening
         */
        std::optional<double> flattening = 0;

        /**
         * @brief Outer radius of the external domain when enabled.
         * @section toml
         * - [main].r_infinity
         */
        std::optional<double> r_infinity = 6.0;

        /**
         * @brief Radius inside which transformations are skipped to avoid singularities.
         * @section toml
         * - [main].r_instability
         */
        std::optional<double> r_instability = 1e-14;
        /**
         * @brief Controls the smoothness/steepness of the core-to-envelope transition.
         * @section toml
         * - [main].core_steepness
         */
        std::optional<double> core_steepness = 1.0;

        /**
         * @brief Continuity order for the core-envelope transition (0 = discontinuous, 1 = C1, 2 = C2).
         * @section toml
         * - [main].continuity_order
         */
        std::optional<size_t> continuity_order = 2;

        /**
         * @brief Boundary attribute id for stellar surface
         * @section toml
         * - [main].surface_bdr_id
         */
        std::optional<size_t> surface_bdr_id = 1;

        /**
         * @brief Boundary attribute id for infinity in kelvin mapping
         * @section toml
         * - [main].inf_bdr_id
         */
        std::optional<size_t> inf_bdr_id = 2;

        /**
         * @brief Material attribute id for the core region
         * @section toml
         * - [main].core_id
         */
        std::optional<size_t> core_id = 1;

        /**
         * @brief Material attribute id for the envelope region
         * @section toml
         * - [main].envelope_id
         */
        std::optional<size_t> envelope_id = 2;

        /**
         * @brief Material attribute id for the external domain (if enabled)
         * @section toml
         * - [main].vacuum_id
         */
        std::optional<size_t> vacuum_id = 3;

        std::optional<OptimizationMethods> optimization_methods = OptimizationMethods{true, true};

        /**
         * @brief Core mapping strategy: legacy "spherified" or conditioned "multi_block".
         *
         * spherified generates a either two or three inscribed cubes then projects them into spheres.
         * multi_block generates a multi-block topology with a single core block and six envelope blocks, then projects the core block into a sphere and the envelope blocks into a spheroid.
         *
         * multi_block is strongly preferred for its ~1000x improved condition number, Spherified is only provided for legacy compatibility.
         *
         * @section toml
         * - [main].core_mapping
         */
        std::optional<std::string> core_mapping = "multi_block";

    };

    inline std::string to_string(const MeshConfig &mesh_config) {
        auto opt_2_string = [](const OptimizationMethods& opt) {
            std::stringstream ss;
            ss << "<OptimizationMethods:";
            if (*opt.tmop) {
                ss << " tmop";
            }
            if (*opt.smoothstep) {
                ss << " smoothstep";
            }
            ss << ">";
            return ss.str();
        };

        std::stringstream ss;

        OptimizationMethods opt = mesh_config.optimization_methods.value_or(OptimizationMethods{false, true});
        std::string opt_string = opt_2_string(opt);

        ss << "MeshConfig:\n";
        ss << std::format("  refinement_levels: {}\n", mesh_config.refinement_levels.value_or(4));
        ss << std::format("  order: {}\n", mesh_config.order.value_or(3));
        ss << std::format("  include_external_domain: {}\n", mesh_config.include_external_domain.value_or(true));
        ss << std::format("  r_core: {}\n", mesh_config.r_core.value_or(0.25));
        ss << std::format("  r_star: {}\n", mesh_config.r_star.value_or(1.0));
        ss << std::format("  flattening: {}\n", mesh_config.flattening.value_or(0.0));
        ss << std::format("  r_infinity: {}\n", mesh_config.r_infinity.value_or(6.0));
        ss << std::format("  r_instability: {}\n", mesh_config.r_instability.value_or(1e-14));
        ss << std::format("  core_steepness: {}\n", mesh_config.core_steepness.value_or(1.0));
        ss << std::format("  continuity_order: {}\n", mesh_config.continuity_order.value_or(2));
        ss << std::format("  surface_bdr_id: {}\n", mesh_config.surface_bdr_id.value_or(1));
        ss << std::format("  inf_bdr_id: {}\n", mesh_config.inf_bdr_id.value_or(2));
        ss << std::format("  core_id: {}\n", mesh_config.core_id.value_or(1));
        ss << std::format("  envelope_id: {}\n", mesh_config.envelope_id.value_or(2));
        ss << std::format("  vacuum_id: {}\n", mesh_config.vacuum_id.value_or(3));
        ss << std::format("  optimization_methods: {}\n", opt_string);
        ss << std::format("  core_mapping: {}\n", mesh_config.core_mapping.value_or("spherified"));

        return ss.str();

    }
}
