#include "stroid/stroid.h"
#include "CLI/CLI.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {
    constexpr double kPi = 3.14159265358979323846;

    class ScopedOutputRedirect {
        std::streambuf* original;
    public:
        ScopedOutputRedirect() : original(std::cout.rdbuf(std::cerr.rdbuf())) {}
        ~ScopedOutputRedirect() { std::cout.rdbuf(original); }
    };

    struct GeometryResults {
        int stellar_elements = 0;
        int vacuum_elements = 0;
        int hanging_faces = 0;
        size_t samples = 0;
        size_t nonpositive_samples = 0;
        double min_signed_det = std::numeric_limits<double>::infinity();
        double max_condition = 0.0;
        double max_face_mismatch = 0.0;
        double max_coordinate_mismatch = 0.0;
        double coordinate_constraint_residual = 0.0;
        double stellar_volume = 0.0;
        double stellar_volume_error = 0.0;
        double surface_radius_error = 0.0;
        double outer_radius_error = 0.0;
    };

    GeometryResults InspectGeometry(stroid::StroidMesh& generated, int grid_points) {
        GeometryResults result;
        auto& mesh = *generated.mesh;
        const int vacuum = static_cast<int>(generated.config.vacuum_id.value());
        const auto& coordinate = *generated.exterior_coordinate->values;
        mfem::Vector independent;
        coordinate.GetTrueDofs(independent);
        mfem::GridFunction reconstructed(generated.exterior_coordinate->space.get());
        reconstructed.SetFromTrueDofs(independent);
        reconstructed -= coordinate;
        result.coordinate_constraint_residual = reconstructed.Normlinf();

        for (int element = 0; element < mesh.GetNE(); ++element) {
            const bool is_vacuum = mesh.GetAttribute(element) == vacuum;
            is_vacuum ? ++result.vacuum_elements : ++result.stellar_elements;
            auto* transformation = mesh.GetElementTransformation(element);
            auto inspect = [&](const mfem::IntegrationPoint& point) {
                transformation->SetIntPoint(&point);
                const auto& jacobian = transformation->Jacobian();
                const double determinant = jacobian.Det();
                const double smallest = jacobian.CalcSingularvalue(2);
                const double largest = jacobian.CalcSingularvalue(0);
                ++result.samples;
                if (!(determinant > 0.0) || !std::isfinite(determinant)) ++result.nonpositive_samples;
                result.min_signed_det = std::min(result.min_signed_det, determinant);
                result.max_condition = std::max(result.max_condition,
                    smallest > 0.0 ? largest / smallest : std::numeric_limits<double>::infinity());
            };
            const auto& quadrature = mfem::IntRules.Get(mfem::Geometry::CUBE, 3 * transformation->Order() + 3);
            for (int q = 0; q < quadrature.GetNPoints(); ++q) {
                const auto& point = quadrature.IntPoint(q);
                inspect(point);
                if (!is_vacuum) result.stellar_volume += point.weight * transformation->Jacobian().Det();
            }
            for (int i = 0; i < grid_points; ++i) {
                for (int j = 0; j < grid_points; ++j) {
                    for (int k = 0; k < grid_points; ++k) {
                        mfem::IntegrationPoint point;
                        point.Set3(static_cast<double>(i) / (grid_points - 1),
                                   static_cast<double>(j) / (grid_points - 1),
                                   static_cast<double>(k) / (grid_points - 1));
                        inspect(point);
                    }
                }
            }
        }
        const double stellar_radius = generated.config.r_star.value();
        const double flattening = generated.config.flattening.value();
        const double analytic_volume = 4.0 * kPi * std::pow(stellar_radius, 3) * (1.0 - flattening) / 3.0;
        result.stellar_volume_error = std::abs(result.stellar_volume - analytic_volume);

        mfem::Vector first(3), second(3);
        for (int face = 0; face < mesh.GetNumFaces(); ++face) {
            const auto information = mesh.GetFaceInformation(face);
            if (!information.IsLocal()) continue;
            if (information.IsNonconformingFine()) ++result.hanging_faces;
            auto* transformation = mesh.GetFaceElementTransformations(face);
            const bool stellar_interface = (transformation->Elem1->Attribute == vacuum) !=
                                           (transformation->Elem2->Attribute == vacuum);
            if (stellar_interface && !information.IsConforming()) {
                throw std::runtime_error("The stellar-vacuum interface is nonconforming.");
            }
            for (int i = 0; i < grid_points; ++i) {
                for (int j = 0; j < grid_points; ++j) {
                    mfem::IntegrationPoint point;
                    point.Set2(static_cast<double>(i) / (grid_points - 1),
                               static_cast<double>(j) / (grid_points - 1));
                    transformation->SetAllIntPoints(&point);
                    const auto& first_point = transformation->Elem1->GetIntPoint();
                    const auto& second_point = transformation->Elem2->GetIntPoint();
                    transformation->Elem1->Transform(first_point, first);
                    transformation->Elem2->Transform(second_point, second);
                    first -= second;
                    result.max_face_mismatch = std::max(result.max_face_mismatch, first.Norml2());
                    const double difference = coordinate.GetValue(transformation->Elem1No, first_point) -
                                              coordinate.GetValue(transformation->Elem2No, second_point);
                    result.max_coordinate_mismatch = std::max(result.max_coordinate_mismatch, std::abs(difference));
                }
            }
        }

        const auto& surface_quadrature = mfem::IntRules.Get(mfem::Geometry::SQUARE, 10);
        for (int boundary = 0; boundary < mesh.GetNBE(); ++boundary) {
            const bool surface = mesh.GetBdrAttribute(boundary) == static_cast<int>(generated.config.surface_bdr_id.value());
            const bool outer = mesh.GetBdrAttribute(boundary) == static_cast<int>(generated.config.inf_bdr_id.value());
            if (!surface && !outer) continue;
            const double radius = surface ? stellar_radius : generated.config.r_infinity.value();
            auto* transformation = mesh.GetBdrElementTransformation(boundary);
            for (int q = 0; q < surface_quadrature.GetNPoints(); ++q) {
                transformation->Transform(surface_quadrature.IntPoint(q), first);
                first(2) /= 1.0 - flattening;
                double& error = surface ? result.surface_radius_error : result.outer_radius_error;
                error = std::max(error, std::abs(first.Norml2() - radius));
            }
        }
        return result;
    }

    struct SolveResults {
        int dofs = 0;
        int true_dofs = 0;
        int iterations = 0;
        bool converged = false;
        double relative_residual = 0.0;
        double l2_error = 0.0;
        double stellar_l2_error = 0.0;
        double vacuum_l2_error = 0.0;
        double h1_error = 0.0;
    };

    SolveResults SolveManufacturedProblem(stroid::StroidMesh& generated, int solution_order) {
        auto& mesh = *generated.mesh;
        const double outer_sixth_power = std::pow(generated.config.r_infinity.value(), 6);
        mfem::FunctionCoefficient exact([outer_sixth_power](const mfem::Vector& point) {
            const double radius_squared = point * point;
            return std::exp(-radius_squared) + 0.1 * std::pow(radius_squared, 3) / outer_sixth_power;
        });
        mfem::VectorFunctionCoefficient gradient(3, [outer_sixth_power](const mfem::Vector& point, mfem::Vector& value) {
            const double radius_squared = point * point;
            value = point;
            value *= -2.0 * std::exp(-radius_squared) + 0.6 * radius_squared * radius_squared / outer_sixth_power;
        });
        mfem::FunctionCoefficient forcing([outer_sixth_power](const mfem::Vector& point) {
            const double radius_squared = point * point;
            return (6.0 - 4.0 * radius_squared) * std::exp(-radius_squared)
                 - 4.2 * radius_squared * radius_squared / outer_sixth_power;
        });
        mfem::H1_FECollection collection(solution_order, 3);
        mfem::FiniteElementSpace space(&mesh, &collection);
        SolveResults result;
        result.dofs = space.GetVSize();
        result.true_dofs = space.GetTrueVSize();
        mfem::Array<int> boundary(mesh.bdr_attributes.Max());
        boundary = 0;
        boundary[static_cast<int>(generated.config.inf_bdr_id.value()) - 1] = 1;
        mfem::Array<int> essential;
        space.GetEssentialTrueDofs(boundary, essential);
        mfem::GridFunction solution(&space);
        solution = 0.0;
        solution.ProjectBdrCoefficient(exact, boundary);
        const int quadrature_order = 2 * solution_order + 3 * generated.config.order.value() + 4;
        const auto& quadrature = mfem::IntRules.Get(mfem::Geometry::CUBE, quadrature_order);
        mfem::LinearForm rhs(&space);
        auto* load = new mfem::DomainLFIntegrator(forcing);
        load->SetIntRule(&quadrature);
        rhs.AddDomainIntegrator(load);
        rhs.Assemble();
        mfem::ConstantCoefficient one(1.0);
        mfem::BilinearForm form(&space);
        auto* diffusion = new mfem::DiffusionIntegrator(one);
        diffusion->SetIntRule(&quadrature);
        form.AddDomainIntegrator(diffusion);
        form.Assemble();
        mfem::OperatorPtr system;
        mfem::Vector independent, system_rhs;
        form.FormLinearSystem(essential, solution, rhs, system, independent, system_rhs);
        mfem::GSSmoother preconditioner(static_cast<mfem::SparseMatrix&>(*system));
        mfem::CGSolver solver;
        solver.SetOperator(*system);
        solver.SetPreconditioner(preconditioner);
        solver.SetRelTol(1.0e-11);
        solver.SetAbsTol(1.0e-14);
        solver.SetMaxIter(2000);
        solver.SetPrintLevel(-1);
        solver.Mult(system_rhs, independent);
        result.converged = solver.GetConverged();
        result.iterations = solver.GetNumIterations();
        mfem::Vector residual(system_rhs.Size());
        system->Mult(independent, residual);
        residual -= system_rhs;
        result.relative_residual = residual.Norml2() / system_rhs.Norml2();
        form.RecoverFEMSolution(independent, rhs, solution);
        const mfem::IntegrationRule* rules[mfem::Geometry::NumGeom]{};
        rules[mfem::Geometry::CUBE] = &quadrature;
        result.l2_error = solution.ComputeL2Error(exact, rules);
        result.h1_error = solution.ComputeH1Error(&exact, &gradient, rules);
        for (int element = 0; element < mesh.GetNE(); ++element) {
            auto* transformation = mesh.GetElementTransformation(element);
            double element_error = 0.0;
            for (int q = 0; q < quadrature.GetNPoints(); ++q) {
                const auto& point = quadrature.IntPoint(q);
                transformation->SetIntPoint(&point);
                const double difference = solution.GetValue(element, point) - exact.Eval(*transformation, point);
                element_error += point.weight * transformation->Weight() * difference * difference;
            }
            if (mesh.GetAttribute(element) == static_cast<int>(generated.config.vacuum_id.value())) {
                result.vacuum_l2_error += element_error;
            } else {
                result.stellar_l2_error += element_error;
            }
        }
        result.stellar_l2_error = std::sqrt(result.stellar_l2_error);
        result.vacuum_l2_error = std::sqrt(result.vacuum_l2_error);
        return result;
    }
}

int main(int argc, char* argv[]) {
    std::vector<int> orders{1, 2};
    std::vector<int> refinements{2, 3};
    std::vector<int> bulk_levels{0};
    int outer_level = -1;
    int solution_order = 1;
    int grid_points = 3;
    double flattening = 0.0;
    double infinity_radius = 6.0;
    std::string output_path;
    CLI::App app{"Compare uniform and graded vacuum meshes, geometry, and an H1 manufactured Poisson solution; TMOP is disabled."};
    app.add_option("--orders", orders, "Geometry orders, comma separated")->delimiter(',')->check(CLI::Range(1, 6));
    app.add_option("--refinements", refinements, "Stellar refinement levels, comma separated")->delimiter(',')->check(CLI::Range(0, 4));
    app.add_option("--bulk-levels", bulk_levels, "Vacuum background target levels, comma separated")->delimiter(',')->check(CLI::Range(0, 4));
    app.add_option("--outer-level", outer_level, "Outer vacuum target; -1 inherits stellar level")->check(CLI::Range(-1, 4));
    app.add_option("--solution-order", solution_order, "H1 polynomial order for the manufactured solve")->check(CLI::Range(1, 4));
    app.add_option("--grid-points", grid_points, "Closed tensor sample grid per coordinate, in addition to quadrature")->check(CLI::Range(2, 9));
    app.add_option("--flattening", flattening, "Spheroidal flattening")->check(CLI::Range(0.0, 0.9));
    app.add_option("--infinity-radius", infinity_radius, "Finite geometric outer radius; stellar radius is one");
    app.add_option("--output", output_path, "New CSV output file; default stdout");
    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError& error) {
        return app.exit(error);
    }
    try {
        if (!std::isfinite(infinity_radius) || infinity_radius <= 1.0) {
            throw std::invalid_argument("Require a finite infinity-radius greater than one.");
        }
        std::ofstream file;
        if (!output_path.empty()) {
            if (std::filesystem::exists(output_path)) throw std::runtime_error("Refusing to overwrite existing output: " + output_path);
            file.open(output_path);
            if (!file) throw std::runtime_error("Could not open output: " + output_path);
        }
        std::ostream& output = output_path.empty() ? std::cout : file;
        output << std::setprecision(17);
        output << "policy,geometry_order,stellar_level,bulk_target,outer_target,solution_order,flattening,r_infinity,elements,stellar_elements,vacuum_elements,vacuum_min_depth,vacuum_max_depth,geometry_dofs,geometry_true_dofs,solution_dofs,solution_true_dofs,hanging_faces,samples,nonpositive_samples,min_signed_det,max_condition,max_face_mismatch,max_coordinate_mismatch,coordinate_constraint_residual,stellar_volume,stellar_volume_error,surface_radius_error,outer_radius_error,l2_error,stellar_l2_error,vacuum_l2_error,h1_error,solver_converged,solver_iterations,relative_residual\n";
        std::cerr << "Evaluating actual FE geometry using signed Jacobians at quadrature points and a closed grid.\n"
                     "Manufactured problem: -Delta u=f, u=exp(-r^2)+0.1*r^6/R^6, exact Dirichlet data only at the outer boundary.\n"
                     "This finite-domain Poisson diagnostic measures discretization error; it does not implement compactified physics.\n";
        bool verified = true;
        for (const int order : orders) {
            for (const int refinement : refinements) {
                std::vector<int> cases{-1};
                cases.insert(cases.end(), bulk_levels.begin(), bulk_levels.end());
                for (const int bulk : cases) {
                    const bool uniform = bulk < 0;
                    const std::string policy = uniform ? "uniform" : "graded";
                    stroid::config::MeshConfig config;
                    config.order = order;
                    config.refinement_levels = refinement;
                    config.flattening = flattening;
                    config.r_infinity = infinity_radius;
                    config.optimization_methods = stroid::config::OptimizationMethods{false, true};
                    if (!uniform) {
                        config.vacuum_refinement_levels = bulk;
                        if (outer_level >= 0) config.vacuum_outer_refinement_levels = outer_level;
                    }
                    std::cerr << "Inspecting " << policy << ", geometry order " << order << ", stellar level " << refinement
                              << ", bulk " << (uniform ? refinement : bulk) << '\n';
                    stroid::StroidMesh generated;
                    {
                        ScopedOutputRedirect redirect;
                        generated = stroid::GenerateMesh(config);
                    }
                    const auto geometry = InspectGeometry(generated, grid_points);
                    const auto stats = stroid::stats::ComputeMeshStats(generated, stroid::stats::MeshStatFeatures::REFINEMENT);
                    const auto solve = SolveManufacturedProblem(generated, solution_order);
                    verified = verified && geometry.nonpositive_samples == 0 && geometry.max_face_mismatch < 1.0e-10
                        && geometry.max_coordinate_mismatch < 1.0e-11 && geometry.coordinate_constraint_residual < 1.0e-11
                        && solve.converged && solve.relative_residual < 1.0e-9 && std::isfinite(solve.h1_error);
                    const int outer = uniform ? refinement : outer_level < 0 ? refinement : outer_level;
                    output << policy << ',' << order << ',' << refinement << ',' << (uniform ? refinement : bulk) << ',' << outer
                           << ',' << solution_order << ',' << flattening << ',' << infinity_radius << ',' << generated.mesh->GetNE()
                           << ',' << geometry.stellar_elements << ',' << geometry.vacuum_elements
                           << ',' << stats.refinement->vacuum.min_depth << ',' << stats.refinement->vacuum.max_depth
                           << ',' << stats.refinement->geometry_dofs << ',' << stats.refinement->geometry_true_dofs
                           << ',' << solve.dofs << ',' << solve.true_dofs << ',' << geometry.hanging_faces
                           << ',' << geometry.samples << ',' << geometry.nonpositive_samples << ',' << geometry.min_signed_det
                           << ',' << geometry.max_condition << ',' << geometry.max_face_mismatch << ',' << geometry.max_coordinate_mismatch
                           << ',' << geometry.coordinate_constraint_residual << ',' << geometry.stellar_volume << ',' << geometry.stellar_volume_error
                           << ',' << geometry.surface_radius_error << ',' << geometry.outer_radius_error << ',' << solve.l2_error
                           << ',' << solve.stellar_l2_error << ',' << solve.vacuum_l2_error << ',' << solve.h1_error
                           << ',' << solve.converged << ',' << solve.iterations << ',' << solve.relative_residual << '\n';
                    output.flush();
                    if (!output) throw std::runtime_error("Failed to write experiment output.");
                }
            }
        }
        if (!verified) {
            std::cerr << "Numerical verification failed: inspect Jacobian, trace, or solver columns.\n";
            return 1;
        }
    } catch (const std::exception& error) {
        std::cerr << "Vacuum refinement experiment failed: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
