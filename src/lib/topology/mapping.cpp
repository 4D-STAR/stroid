#include "stroid/topology/mapping.h"
#include "stroid/exceptions/exceptions.h"
#include <cmath>
#include <algorithm>
#include <array>
#include <utility>
#include <format>
#include <string>

namespace {
    template<int n, int k>
    consteval int nCr() {
        if constexpr (k > n) {
            return 0;
        } else {
            if constexpr (constexpr int kk = (k * 2 > n) ? (n - k) : k; kk == 0) {
                return 1;
            } else {
                int result = n;

                for (int i = 2; i <= kk; ++i) {
                    result *= (n - i + 1);
                    result /= i;
                }

                return result;
            }
        }
    }

    template <int n>
    double GeneralizedSmoothstep(const double x) {
        if (x <= 0.0) return 0.0;
        if (x >= 1.0) return 1.0;

        double sum = 0.0;

        auto compute_term = [&]<std::size_t k>(std::integral_constant<std::size_t, k>) {
            return nCr<n + k, k>() * std::pow(1.0 - x, k);
        };

        auto unroller = [&]<std::size_t... ks>(std::index_sequence<ks...>) {
            return (compute_term(std::integral_constant<std::size_t, ks>{}) + ...);
        };

        sum = unroller(std::make_index_sequence<n + 1>{});
        return sum * std::pow(x, n + 1);
    }

    template <std::size_t... Is>
    constexpr auto make_smoothstep_dispatch_table(std::index_sequence<Is...>) {
        return std::array<double(*)(double), sizeof...(Is)>{
            &GeneralizedSmoothstep<Is + 1>...
        };
    }

    constexpr int MAX_SMOOTHSTEP_ORDER = 10;
    constexpr auto smoothstep_dispatch = make_smoothstep_dispatch_table(
        std::make_index_sequence<MAX_SMOOTHSTEP_ORDER>{}
    );
}

namespace stroid::topology {
    void ApplyEquiangular(mfem::Vector &pos) {
        const double x = pos(0);
        const double y = pos(1);
        const double z = pos(2);

        const double absX = std::abs(x);
        const double absY = std::abs(y);
        const double absZ = std::abs(z);

        const double maxAbs = std::max({absX, absY, absZ});

        if (maxAbs < 1e-14) return;

        if (absX == maxAbs) {
            pos(1) = x * std::tan(M_PI / 4.0 * (y/x));
            pos(2) = x * std::tan(M_PI / 4.0 * (z/x));
        } else if (absY == maxAbs) {
            pos(0) = y * std::tan(M_PI / 4.0 * (x/y));
            pos(2) = y * std::tan(M_PI / 4.0 * (z/y));
        } else { // absZ == maxAbs
            pos(0) = z * std::tan(M_PI / 4.0 * (x/z));
            pos(1) = z * std::tan(M_PI / 4.0 * (y/z));
        }
    }

    void ApplySpheroidal(mfem::Vector &pos, const fourdst::config::Config<config::MeshConfig> &config) {
        pos(2) *= (1.0 - config->flattening.value());
    }

    void TransformPoint(mfem::Vector &pos, const fourdst::config::Config<config::MeshConfig> &config, int attribute_id) {
        double X = pos(0);
        double Y = pos(1);
        double Z = pos(2);

        double maxAbs = std::max({std::abs(X), std::abs(Y), std::abs(Z)});
        const bool multi_block = config->core_mapping.value_or("spherified") == "multi_block";
        const double inner_radius = config->r_core.value() / 2.0;
        if (multi_block && maxAbs <= inner_radius) {
            pos /= std::sqrt(3.0);
            ApplySpheroidal(pos, config);
            return;
        }
        if (!multi_block && maxAbs < 1e-14) return;

        double cx = X / maxAbs;
        double cy = Y / maxAbs;
        double cz = Z / maxAbs;

        double sx = cx * std::sqrt(1.0 - cy*cy/2.0 - cz*cz/2.0 + cy*cy*cz*cz/3.0);
        double sy = cy * std::sqrt(1.0 - cx*cx/2.0 - cz*cz/2.0 + cx*cx*cz*cz/3.0);
        double sz = cz * std::sqrt(1.0 - cx*cx/2.0 - cy*cy/2.0 + cx*cx*cy*cy/3.0);

        mfem::Vector unit_dir(3);
        unit_dir(0) = sx;
        unit_dir(1) = sy;
        unit_dir(2) = sz;

        if (maxAbs <= config->r_core.value()) {
            if (multi_block) {
                const double t = (maxAbs - inner_radius) / inner_radius;
                const double inner_scale = inner_radius / std::sqrt(3.0);
                pos(0) = (1.0 - t) * inner_scale * cx + t * config->r_core.value() * sx;
                pos(1) = (1.0 - t) * inner_scale * cy + t * config->r_core.value() * sy;
                pos(2) = (1.0 - t) * inner_scale * cz + t * config->r_core.value() * sz;
                ApplySpheroidal(pos, config);
                return;
            }
            double nx = X / config->r_core.value();
            double ny = Y / config->r_core.value();
            double nz = Z / config->r_core.value();

            pos(0) = nx * std::sqrt(1.0 - ny*ny/2.0 - nz*nz/2.0 + ny*ny*nz*nz/3.0);
            pos(1) = ny * std::sqrt(1.0 - nx*nx/2.0 - nz*nz/2.0 + nx*nx*nz*nz/3.0);
            pos(2) = nz * std::sqrt(1.0 - nx*nx/2.0 - ny*ny/2.0 + nx*nx*ny*ny/3.0);

            pos *= config->r_core.value();

            ApplySpheroidal(pos, config);
            return;
        }

        if (maxAbs <= config->r_star.value()) {
            const double xi = (maxAbs - config->r_core.value()) / (config->r_star.value() - config->r_core.value());
            const double r_phys = config->r_core.value() + xi * (config->r_star.value() - config->r_core.value());

            pos = unit_dir;
            pos *= r_phys;

            ApplySpheroidal(pos, config);
        } else {
            pos = unit_dir;
            pos *= maxAbs;

            ApplySpheroidal(pos, config);
        }
    }

    double ComputeExteriorCoordinate(
        const mfem::Vector& logical_position,
        const int attribute,
        const fourdst::config::Config<config::MeshConfig>& config
    ) {
        if (!config->include_external_domain.value() || attribute != static_cast<int>(config->vacuum_id.value())) return 0.0;

        const double logical_radius = std::max({
            std::abs(logical_position(0)),
            std::abs(logical_position(1)),
            std::abs(logical_position(2))
        });

        const double r_star = config->r_star.value();
        const double r_infinity = config->r_infinity.value();
        const double coordinate = (logical_radius - r_star) / (r_infinity - r_star);
        constexpr double tolerance = 64.0 * std::numeric_limits<double>::epsilon();

        if (coordinate < -tolerance || coordinate > 1.0 + tolerance) {
            throw std::runtime_error("Logical exterior coordinate lies outside [0, 1].");
        }

        if (std::abs(coordinate) <= tolerance) return 0.0;
        if (std::abs(coordinate - 1.0) <= tolerance) return 1.0;
        return coordinate;
    }

}
