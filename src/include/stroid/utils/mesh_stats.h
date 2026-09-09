#pragma once

#include "mfem.hpp"
#include "stroid/utils/types.h"
#include "stroid/config/config.h"

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace stroid::stats {
    enum class MeshStatFeatures : uint32_t {
        NONE                 = 0u,
        RADIUS               = 1u << 0,
        AXES                 = 1u << 1,
        ELLIPTICITY          = 1u << 2,
        BOWING               = 1u << 3,
        CONFORMITY           = 1u << 4,
        JACOBIAN             = 1u << 5,
        VOLUME_AREA          = 1u << 6,
        ELEMENT_COUNT        = 1u << 7,
        MESH_SIZE            = 1u << 8,
        OUTER_BOUNDS         = 1u << 9,
        CENTROID             = 1u << 10,
        CONFIG_META          = 1u << 11,
        BOUNDING_BOX         = 1u << 12,
        REFINEMENT           = 1u << 13,
    };

    constexpr MeshStatFeatures operator|(MeshStatFeatures lhs, MeshStatFeatures rhs) {
        return static_cast<MeshStatFeatures>(static_cast<uint32_t>(lhs) | static_cast<uint32_t>(rhs));
    }

    constexpr MeshStatFeatures operator&(MeshStatFeatures lhs, MeshStatFeatures rhs) {
        return static_cast<MeshStatFeatures>(static_cast<uint32_t>(lhs) & static_cast<uint32_t>(rhs));
    }

    constexpr bool has_feature(MeshStatFeatures feature, MeshStatFeatures set) {
        return (static_cast<uint32_t>(set) & static_cast<uint32_t>(feature)) != 0u;
    }

    inline constexpr MeshStatFeatures MESH_STAT_DEFAULT =
        MeshStatFeatures::RADIUS | MeshStatFeatures::AXES | MeshStatFeatures::ELLIPTICITY |
        MeshStatFeatures::CONFORMITY | MeshStatFeatures::CONFIG_META | MeshStatFeatures::REFINEMENT;

    inline constexpr auto MESH_STAT_ALL = static_cast<MeshStatFeatures>(0xFFFFFFFFu);

    struct RadiusStats {
        double min = 0, max = 0, mean = 0, stddev = 0;
        long n_samples = 0;
    };

    struct AxisStats {
        double semi_major = 0;
        double semi_minor = 0;
    };

    struct EllipticityStats {
        double flattening = 0;
        double polar_equatorial = 1;
        double radius_uniformity = 1;
    };

    struct BowingStats {
        double max_inward = 0;
        double max_outward = 0;
        double rms = 0;
        double worst_at_radius = 0;
    };

    struct ConformityStats {
        bool conforming = true;
        bool hierarchy_enabled = false;
        // Fine patches are counted once; their coarse master faces are excluded.
        long n_nonconforming_faces = 0;
    };

    struct RegionRefinementStats {
        // An absent region has both depths set to -1.
        int min_depth = -1;
        int max_depth = -1;
    };

    struct RefinementStats {
        RegionRefinementStats all;
        RegionRefinementStats core;
        RegionRefinementStats envelope;
        RegionRefinementStats vacuum;
        // Scalar nodal counts, independent of the coordinate vector dimension.
        long geometry_dofs = 0;
        long geometry_true_dofs = 0;
    };

    struct JacobianStats {
        double detJ_min;
        double detJ_max;
        long   n_flipped;
        double min_detJ_ratio;
        double worst_ratio_at_radius;
        double detJ_min_at_radius;
        long   n_elements;
    };
    struct VolumeAreaStats {
        double stellar_volume = 0, surface_area = 0;
        double analytic_volume = 0, analytic_area = 0;
    };

    struct ElementCounts {
        long total = 0, core = 0, envelope = 0, vacuum = 0, other = 0;
        long n_vertices = 0;
    };

    struct MeshSizeStats {
        double h_min = 0, h_max = 0, h_mean = 0, h_stddev = 0;
    };

    struct OuterBoundsStats {
        double min = 0, max = 0, mean = 0;
        long n_samples = 0;
    };

    struct CentroidStats {
        double x = 0, y = 0, z = 0, offset = 0;
    };

    struct ConfigMeta {
        double r_core = 0, r_star = 0, flattening = 0, r_infinity = 0;
        int geom_order = 0;
        size_t refinement_levels = 0;
        bool has_external_domain = true;
    };

    struct BoundingBox {
        double xMin = 0, xMax = 0, yMin = 0, yMax = 0, zMin = 0, zMax = 0;
        bool valid = false;

        [[nodiscard]] double dx() const {return xMax - xMin;}
        [[nodiscard]] double dy() const {return yMax - yMin;}
        [[nodiscard]] double dz() const {return zMax - zMin;}
        [[nodiscard]] double diag() const {
            const double a = dx(), b = dy(), c = dz();
            return std::sqrt(a*a + b*b + c*c);
        }
    };

    struct BoundingBoxStats {
        BoundingBox core;
        BoundingBox star;
        BoundingBox vacuum;
    };

    struct MeshStats {
        MeshStatFeatures computed = MeshStatFeatures::NONE;
        std::optional<RadiusStats>                      radius;
        std::optional<AxisStats>                        axes;
        std::optional<EllipticityStats>                 ellipticity;
        std::optional<BowingStats>                      bowing;
        std::optional<ConformityStats>                  conformity;
        std::optional<RefinementStats>                  refinement;
        std::optional<JacobianStats>                    jacobian;
        std::optional<JacobianStats>                    jacobian_stellar;
        std::optional<JacobianStats>                    jacobian_vacuum;
        std::optional<VolumeAreaStats>                  volume;
        std::optional<ElementCounts>                    element_counts;
        std::optional<MeshSizeStats>                    mesh_size;
        std::optional<OuterBoundsStats>                 outer_bounds;
        std::optional<CentroidStats>                    centroid;
        std::optional<ConfigMeta>                       config_meta;
        std::optional<BoundingBoxStats>                 bounding_box;

        std::vector<std::string>                        warnings;
        std::vector<std::string>                        errors;
    };

    MeshStats ComputeMeshStats(const StroidMesh& sm, MeshStatFeatures features = MESH_STAT_DEFAULT, int sample_order = -1);

    std::string to_string(const MeshStats& s);

    inline std::ostream& operator<<(std::ostream& os, const MeshStats& s) {
        return os << to_string(s);
    }
}

template <>
struct std::formatter<stroid::stats::MeshStats, char> {
    static constexpr auto parse(const std::format_parse_context& ctx) {
        return ctx.begin();
    }

    static auto format(const stroid::stats::MeshStats &s, std::format_context& ctx) {
        return std::format_to(ctx.out(), "{}", stroid::stats::to_string(s));
    }
};
