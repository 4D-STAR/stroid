#include "stroid/stroid.h"
#include "CLI/CLI.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace {
    constexpr double infinity = std::numeric_limits<double>::infinity();

    struct SampleLocation {
        int element = -1;
        mfem::IntegrationPoint point;
        std::string source = "none";
    };

    struct ConditioningStats {
        int elements = 0;
        size_t samples = 0;
        size_t nonpositive_samples = 0;
        double min_det = infinity;
        double min_sigma = infinity;
        double max_condition = 0.0;
        double min_scaled_jacobian = infinity;
        double contraction_boundary = infinity;
        SampleLocation det_location;
        SampleLocation condition_location;
        SampleLocation contraction_location;
    };

    class ScopedOutputRedirect {
        std::streambuf* original;
    public:
        ScopedOutputRedirect() : original(std::cout.rdbuf(std::cerr.rdbuf())) {}
        ~ScopedOutputRedirect() { std::cout.rdbuf(original); }
    };

    std::vector<mfem::IntegrationPoint> ClosedSamples(int grid_points) {
        std::vector<mfem::IntegrationPoint> points;
        for (int i = 0; i < grid_points; ++i) {
            for (int j = 0; j < grid_points; ++j) {
                for (int k = 0; k < grid_points; ++k) {
                    mfem::IntegrationPoint point;
                    point.Set3(static_cast<double>(i) / (grid_points - 1),
                               static_cast<double>(j) / (grid_points - 1),
                               static_cast<double>(k) / (grid_points - 1));
                    points.push_back(point);
                }
            }
        }
        for (const double offset : {0.005, 0.010885670927, 0.02}) {
            for (int corner = 0; corner < 8; ++corner) {
                mfem::IntegrationPoint point;
                point.Set3((corner & 1) ? 1.0 - offset : offset,
                           (corner & 2) ? 1.0 - offset : offset,
                           (corner & 4) ? 1.0 - offset : offset);
                points.push_back(point);
            }
        }
        return points;
    }

    double ColumnNorm(const mfem::DenseMatrix& matrix, int column) {
        double norm_squared = 0.0;
        for (int row = 0; row < 3; ++row) {
            norm_squared += matrix(row, column) * matrix(row, column);
        }
        return std::sqrt(norm_squared);
    }

    double ContractionBoundary(const mfem::DenseMatrix& jacobian,
                               const mfem::DenseMatrix& direction) {
        std::array<double, 4> coefficients{};
        mfem::DenseMatrix mixed(3);
        for (int mask = 0; mask < 8; ++mask) {
            int degree = 0;
            for (int column = 0; column < 3; ++column) {
                const bool use_direction = (mask & (1 << column)) != 0;
                degree += use_direction;
                for (int row = 0; row < 3; ++row) {
                    mixed(row, column) = use_direction ? direction(row, column) : jacobian(row, column);
                }
            }
            coefficients[degree] += mixed.Det();
        }
        if (!(coefficients[0] > 0.0)) return 0.0;

        auto polynomial = [&](double x) {
            return ((coefficients[3] * x + coefficients[2]) * x + coefficients[1]) * x + coefficients[0];
        };
        std::vector<double> breaks{0.0, 1.0};
        auto add_break = [&](double x) {
            if (std::isfinite(x) && x > 0.0 && x < 1.0) breaks.push_back(x);
        };
        const double a = 3.0 * coefficients[3];
        const double b = 2.0 * coefficients[2];
        const double c = coefficients[1];
        if (a == 0.0) {
            if (b != 0.0) add_break(-c / b);
        } else {
            const double discriminant = b * b - 4.0 * a * c;
            if (discriminant >= 0.0) {
                const double q = -0.5 * (b + std::copysign(std::sqrt(discriminant), b));
                if (q == 0.0) {
                    add_break(-b / (2.0 * a));
                } else {
                    add_break(q / a);
                    add_break(c / q);
                }
            }
        }
        std::sort(breaks.begin(), breaks.end());
        const double scale = std::abs(coefficients[0]) + std::abs(coefficients[1])
                           + std::abs(coefficients[2]) + std::abs(coefficients[3]);
        const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() * scale;
        for (size_t i = 1; i < breaks.size(); ++i) {
            double right = breaks[i];
            const double value = polynomial(right);
            if (value > tolerance) continue;
            if (std::abs(value) <= tolerance) return right;
            double left = breaks[i - 1];
            for (int iteration = 0; iteration < 64; ++iteration) {
                const double middle = 0.5 * (left + right);
                if (polynomial(middle) > 0.0) left = middle;
                else right = middle;
            }
            return right;
        }
        return infinity;
    }

    std::unique_ptr<mfem::GridFunction> BuildContractionProbe(
        stroid::StroidMesh& mesh, mfem::FiniteElementSpace& space, double stellar_radius
    ) {
        auto values = std::make_unique<mfem::GridFunction>(&space);
        *values = 0.0;
        std::vector<bool> processed(space.GetNDofs(), false);
        mfem::Array<int> dofs;
        mfem::Vector physical(3), logical(3);
        for (int element = 0; element < mesh.mesh->GetNE(); ++element) {
            const auto& nodes = space.GetFE(element)->GetNodes();
            space.GetElementDofs(element, dofs);
            auto* physical_transform = mesh.mesh->GetElementTransformation(element);
            auto* logical_transform = mesh.reference_mesh->GetElementTransformation(element);
            for (int local = 0; local < dofs.Size(); ++local) {
                const int dof = dofs[local] >= 0 ? dofs[local] : -1 - dofs[local];
                if (processed[dof]) continue;
                physical_transform->Transform(nodes.IntPoint(local), physical);
                logical_transform->Transform(nodes.IntPoint(local), logical);
                const double radius = physical.Norml2();
                const double logical_radius = std::max({std::abs(logical(0)), std::abs(logical(1)), std::abs(logical(2))});
                const double fraction = std::min(logical_radius / stellar_radius, 1.0);
                for (int component = 0; component < 3; ++component) {
                    (*values)(space.DofToVDof(dof, component)) = radius > 0.0
                        ? -stellar_radius * fraction * fraction * physical(component) / radius : 0.0;
                }
                processed[dof] = true;
            }
        }
        return values;
    }

    std::map<int, ConditioningStats> InspectMesh(stroid::StroidMesh& mesh, int order,
        int grid_points, bool contraction_probe, int probe_order) {
        std::map<int, ConditioningStats> result;
        const auto closed_samples = ClosedSamples(grid_points);
        mfem::H1_FECollection probe_collection(probe_order, 3);
        mfem::FiniteElementSpace probe_space(mesh.mesh.get(), &probe_collection, 3);
        std::unique_ptr<mfem::GridFunction> probe;
        if (contraction_probe) probe = BuildContractionProbe(mesh, probe_space, mesh.config.r_star.value());
        mfem::DenseMatrix probe_values, probe_shape, direction(3);
        mfem::Array<int> probe_dofs;

        for (int element = 0; element < mesh.mesh->GetNE(); ++element) {
            const int attribute = mesh.mesh->GetAttribute(element);
            ++result[attribute].elements;
            ++result[0].elements;
            auto* transform = mesh.mesh->GetElementTransformation(element);
            const bool inspect_probe = contraction_probe && attribute != static_cast<int>(mesh.config.vacuum_id.value());
            if (inspect_probe) {
                probe_space.GetElementDofs(element, probe_dofs);
                probe_values.SetSize(3, probe_dofs.Size());
                probe_shape.SetSize(probe_dofs.Size(), 3);
                for (int local = 0; local < probe_dofs.Size(); ++local) {
                    const int dof = probe_dofs[local] >= 0 ? probe_dofs[local] : -1 - probe_dofs[local];
                    for (int component = 0; component < 3; ++component) {
                        probe_values(component, local) = (*probe)(probe_space.DofToVDof(dof, component));
                    }
                }
            }
            auto inspect_point = [&](const mfem::IntegrationPoint& point, const std::string& source) {
                transform->SetIntPoint(&point);
                const mfem::DenseMatrix& jacobian = transform->Jacobian();
                const double determinant = jacobian.Det();
                const double sigma_min = jacobian.CalcSingularvalue(2);
                const double sigma_max = jacobian.CalcSingularvalue(0);
                const double condition = sigma_min > 0.0 ? sigma_max / sigma_min : infinity;
                const double denominator = ColumnNorm(jacobian, 0) * ColumnNorm(jacobian, 1) * ColumnNorm(jacobian, 2);
                const double scaled_jacobian = denominator > 0.0 ? determinant / denominator : 0.0;
                double boundary = infinity;
                if (inspect_probe) {
                    probe_space.GetFE(element)->CalcDShape(point, probe_shape);
                    mfem::Mult(probe_values, probe_shape, direction);
                    boundary = ContractionBoundary(jacobian, direction);
                }
                for (const int region : {0, attribute}) {
                    auto& stats = result[region];
                    ++stats.samples;
                    if (!(determinant > 0.0)) ++stats.nonpositive_samples;
                    if (determinant < stats.min_det) {
                        stats.min_det = determinant;
                        stats.det_location = {element, point, source};
                    }
                    stats.min_sigma = std::min(stats.min_sigma, sigma_min);
                    stats.min_scaled_jacobian = std::min(stats.min_scaled_jacobian, scaled_jacobian);
                    if (condition > stats.max_condition) {
                        stats.max_condition = condition;
                        stats.condition_location = {element, point, source};
                    }
                    if (boundary < stats.contraction_boundary) {
                        stats.contraction_boundary = boundary;
                        stats.contraction_location = {element, point, source};
                    }
                }
            };
            const auto& quadrature = mfem::IntRules.Get(transform->GetGeometryType(), 2 * order + 4);
            for (int point = 0; point < quadrature.GetNPoints(); ++point) {
                inspect_point(quadrature.IntPoint(point), "quadrature");
            }
            for (const auto& point : closed_samples) inspect_point(point, "closed_grid_and_corner_probes");
        }
        return result;
    }

    void WriteLocation(std::ostream& output, const SampleLocation& location, stroid::StroidMesh& mesh) {
        output << ',' << location.element << ',' << location.source;
        if (location.element < 0) {
            output << ",nan,nan,nan,nan,nan,nan,nan,nan,nan";
            return;
        }
        mfem::Vector physical(3), logical(3);
        mesh.mesh->GetElementTransformation(location.element)->Transform(location.point, physical);
        mesh.reference_mesh->GetElementTransformation(location.element)->Transform(location.point, logical);
        output << ',' << location.point.x << ',' << location.point.y << ',' << location.point.z;
        for (int component = 0; component < 3; ++component) output << ',' << physical(component);
        for (int component = 0; component < 3; ++component) output << ',' << logical(component);
    }
}

int main(int argc, char** argv) {
    std::vector<int> orders{1, 2, 3, 4, 5, 6};
    std::vector<int> refinements{0, 1, 2};
    std::vector<std::string> mappings{"spherified", "multi_block"};
    std::string output_path;
    int grid_points = 5;
    int probe_order = 3;
    double core_radius = 0.25;
    double infinity_radius = 5.0;
    double flattening = 0.0;
    bool no_external = false;
    bool contraction_probe = false;
    CLI::App app{"Compare signed Jacobians and conditioning of the actual high-order STROID mesh; TMOP is disabled."};
    app.add_option("--orders", orders, "Geometry orders, comma separated")->delimiter(',')->check(CLI::Range(1, 8));
    app.add_option("--refinements", refinements, "Uniform refinement levels, comma separated")->delimiter(',')->check(CLI::Range(0, 3));
    app.add_option("--mappings", mappings, "Core mappings, comma separated")->delimiter(',')->check(CLI::IsMember({"spherified", "multi_block"}));
    app.add_option("--grid-points", grid_points, "Closed tensor grid points per coordinate, plus near-corner probes")->check(CLI::Range(2, 15));
    app.add_option("--core-radius", core_radius, "Core radius; stellar radius is one");
    app.add_option("--infinity-radius", infinity_radius, "Outer reference radius");
    app.add_option("--flattening", flattening, "Spheroidal flattening");
    app.add_option("--output", output_path, "New CSV output file; defaults to stdout");
    app.add_flag("--no-external", no_external, "Omit exterior domain");
    app.add_flag("--contraction-probe", contraction_probe, "Inspect an interpolated unit logical-radius-squared radial contraction in stellar elements");
    app.add_option("--probe-order", probe_order, "H1 displacement order for the optional contraction probe")->check(CLI::Range(1, 8));
    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError& error) {
        return app.exit(error);
    }

    try {
        if (!std::isfinite(core_radius) || core_radius <= 0.0 || core_radius >= 1.0
            || !std::isfinite(infinity_radius) || infinity_radius <= 1.0
            || !std::isfinite(flattening) || flattening < 0.0 || flattening >= 1.0) {
            throw std::invalid_argument("Require 0 < core-radius < 1 < infinity-radius and 0 <= flattening < 1.");
        }
        std::ofstream file;
        if (!output_path.empty()) {
            if (std::filesystem::exists(output_path)) throw std::runtime_error("Refusing to overwrite existing output: " + output_path);
            file.open(output_path);
            if (!file) throw std::runtime_error("Could not open output: " + output_path);
        }
        std::ostream& output = output_path.empty() ? std::cout : file;
        output << std::setprecision(17);
        output << "mapping,order,refinement,r_core,r_star,r_infinity,flattening,external,grid_points,quadrature_order,probe_order,attribute,elements,samples,nonpositive_samples,min_signed_det,min_sigma,max_condition,min_scaled_jacobian,contraction_boundary_up_to_one";
        for (const std::string prefix : {"det", "condition", "contraction"}) {
            output << ',' << prefix << "_element," << prefix << "_source," << prefix << "_xi," << prefix << "_eta," << prefix << "_zeta," << prefix << "_x," << prefix << "_y," << prefix << "_z," << prefix << "_logical_x," << prefix << "_logical_y," << prefix << "_logical_z";
        }
        output << '\n';
        std::cerr << "Sampling actual FE geometry, not the analytical map. Attribute 0 aggregates all regions.\n"
                     "Signed determinants and scaled Jacobians retain orientation; inf boundary means no sampled root through alpha=1.\n"
                     "The optional contraction probe is a diagnostic field, not a Newton correction or a production exterior extension.\n";
        for (const auto& mapping : mappings) {
            for (const int order : orders) {
                for (const int refinement : refinements) {
                    stroid::config::MeshConfig config;
                    config.core_mapping = mapping;
                    config.order = order;
                    config.refinement_levels = refinement;
                    config.r_core = core_radius;
                    config.r_star = 1.0;
                    config.r_infinity = infinity_radius;
                    config.flattening = flattening;
                    config.include_external_domain = !no_external;
                    config.optimization_methods = stroid::config::OptimizationMethods{false, true};
                    std::cerr << "Inspecting " << mapping << ", order " << order << ", refinement " << refinement << '\n';
                    stroid::StroidMesh mesh;
                    {
                        ScopedOutputRedirect redirect;
                        mesh = stroid::GenerateMesh(config);
                    }
                    const auto regions = InspectMesh(mesh, order, grid_points, contraction_probe, probe_order);
                    for (const auto& [attribute, stats] : regions) {
                        output << mapping << ',' << order << ',' << refinement << ',' << core_radius << ",1," << infinity_radius << ',' << flattening << ',' << !no_external << ',' << grid_points << ',' << 2 * order + 4 << ',' << (contraction_probe ? probe_order : 0) << ',' << attribute << ',' << stats.elements << ',' << stats.samples << ',' << stats.nonpositive_samples << ',' << stats.min_det << ',' << stats.min_sigma << ',' << stats.max_condition << ',' << stats.min_scaled_jacobian << ',' << stats.contraction_boundary;
                        WriteLocation(output, stats.det_location, mesh);
                        WriteLocation(output, stats.condition_location, mesh);
                        WriteLocation(output, stats.contraction_location, mesh);
                        output << '\n';
                    }
                    output.flush();
                    if (!output) throw std::runtime_error("Failed to write experiment output.");
                }
            }
        }
    } catch (const std::exception& error) {
        std::cerr << "Geometry quality experiment failed: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
