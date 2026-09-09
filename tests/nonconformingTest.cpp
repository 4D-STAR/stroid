#include <gtest/gtest.h>

#include "stroid/stroid.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <map>
#include <string>

namespace {
    constexpr double kPi = 3.14159265358979323846;

    stroid::config::MeshConfig Configuration(int order = 2, int stellar_level = 2,
                                             int bulk_level = 0, int outer_level = 0) {
        stroid::config::MeshConfig config;
        config.order = order;
        config.refinement_levels = stellar_level;
        config.vacuum_refinement_levels = bulk_level;
        config.vacuum_outer_refinement_levels = outer_level;
        config.optimization_methods = stroid::config::OptimizationMethods{false, true};
        return config;
    }

    std::map<int, int> CountAttributes(const mfem::Mesh& mesh, bool boundary = false) {
        std::map<int, int> counts;
        for (int element = 0; element < (boundary ? mesh.GetNBE() : mesh.GetNE()); ++element) {
            ++counts[boundary ? mesh.GetBdrAttribute(element) : mesh.GetAttribute(element)];
        }
        return counts;
    }

    void ExpectConstrainedField(mfem::GridFunction& values) {
        mfem::Vector independent;
        values.GetTrueDofs(independent);
        mfem::GridFunction reconstructed(values.FESpace());
        reconstructed.SetFromTrueDofs(independent);
        reconstructed -= values;
        EXPECT_LT(reconstructed.Normlinf(), 5.0e-12);
    }

    void ExpectGeometryAndCoordinateTraces(stroid::StroidMesh& generated, bool require_hanging_faces = true) {
        mfem::Mesh& mesh = *generated.mesh;
        ASSERT_NE(generated.exterior_coordinate, nullptr);
        mfem::GridFunction& coordinate = *generated.exterior_coordinate->values;
        ASSERT_EQ(coordinate.FESpace()->GetMesh(), &mesh);
        ExpectConstrainedField(*mesh.GetNodes());
        ExpectConstrainedField(coordinate);
        const int vacuum = static_cast<int>(generated.config.vacuum_id.value());
        int hanging_faces = 0;
        int stellar_faces = 0;
        double geometry_error = 0.0;
        double coordinate_error = 0.0;
        double stellar_trace_error = 0.0;
        mfem::Vector first(3), second(3);

        for (int face = 0; face < mesh.GetNumFaces(); ++face) {
            const auto info = mesh.GetFaceInformation(face);
            if (!info.IsLocal()) continue;
            auto* transformation = mesh.GetFaceElementTransformations(face);
            ASSERT_NE(transformation->Elem1, nullptr);
            ASSERT_NE(transformation->Elem2, nullptr);
            const bool first_vacuum = transformation->Elem1->Attribute == vacuum;
            const bool second_vacuum = transformation->Elem2->Attribute == vacuum;
            const bool stellar_interface = first_vacuum != second_vacuum;
            if (info.IsNonconformingFine()) {
                ++hanging_faces;
                EXPECT_EQ(first_vacuum, second_vacuum);
            }
            if (stellar_interface) {
                ++stellar_faces;
                EXPECT_TRUE(info.IsConforming()) << "Stellar interface face " << face;
            }
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    mfem::IntegrationPoint point;
                    point.Set2(i / 3.0, j / 3.0);
                    transformation->SetAllIntPoints(&point);
                    const auto& first_point = transformation->Elem1->GetIntPoint();
                    const auto& second_point = transformation->Elem2->GetIntPoint();
                    transformation->Elem1->Transform(first_point, first);
                    transformation->Elem2->Transform(second_point, second);
                    first -= second;
                    geometry_error = std::max(geometry_error, first.Norml2());
                    const double first_value = coordinate.GetValue(transformation->Elem1No, first_point);
                    const double second_value = coordinate.GetValue(transformation->Elem2No, second_point);
                    coordinate_error = std::max(coordinate_error, std::abs(first_value - second_value));
                    if (stellar_interface) {
                        stellar_trace_error = std::max({stellar_trace_error, std::abs(first_value), std::abs(second_value)});
                    }
                }
            }
        }
        if (require_hanging_faces) EXPECT_GT(hanging_faces, 0);
        EXPECT_GT(stellar_faces, 0);
        EXPECT_LT(geometry_error, 5.0e-12);
        EXPECT_LT(coordinate_error, 5.0e-12);
        EXPECT_LT(stellar_trace_error, 5.0e-12);

        int outer_faces = 0;
        double outer_trace_error = 0.0;
        for (int boundary = 0; boundary < mesh.GetNBE(); ++boundary) {
            if (mesh.GetBdrAttribute(boundary) != static_cast<int>(generated.config.inf_bdr_id.value())) continue;
            ++outer_faces;
            auto* transformation = mesh.GetBdrFaceTransformations(boundary);
            const auto& quadrature = mfem::IntRules.Get(mfem::Geometry::SQUARE, 6);
            for (int q = 0; q < quadrature.GetNPoints(); ++q) {
                transformation->SetAllIntPoints(&quadrature.IntPoint(q));
                outer_trace_error = std::max(outer_trace_error,
                    std::abs(coordinate.GetValue(transformation->Elem1No, transformation->Elem1->GetIntPoint()) - 1.0));
            }
        }
        EXPECT_GT(outer_faces, 0);
        EXPECT_LT(outer_trace_error, 5.0e-12);

        double minimum = 1.0;
        double maximum = 0.0;
        double interior_error = 0.0;
        for (int element = 0; element < mesh.GetNE(); ++element) {
            const auto& quadrature = mfem::IntRules.Get(mfem::Geometry::CUBE, 6);
            for (int q = 0; q < quadrature.GetNPoints(); ++q) {
                const double value = coordinate.GetValue(element, quadrature.IntPoint(q));
                ASSERT_TRUE(std::isfinite(value));
                if (mesh.GetAttribute(element) == vacuum) {
                    minimum = std::min(minimum, value);
                    maximum = std::max(maximum, value);
                } else {
                    interior_error = std::max(interior_error, std::abs(value));
                }
            }
        }
        EXPECT_GE(minimum, -5.0e-12);
        EXPECT_LE(maximum, 1.0 + 5.0e-12);
        EXPECT_LT(interior_error, 5.0e-12);
    }

    void ExpectPositiveJacobians(mfem::Mesh& mesh, int excluded_attribute = -1) {
        double minimum = std::numeric_limits<double>::infinity();
        int minimum_element = -1;
        for (int element = 0; element < mesh.GetNE(); ++element) {
            if (mesh.GetAttribute(element) == excluded_attribute) continue;
            auto* transformation = mesh.GetElementTransformation(element);
            auto inspect = [&](const mfem::IntegrationPoint& point) {
                transformation->SetIntPoint(&point);
                const double determinant = transformation->Jacobian().Det();
                ASSERT_TRUE(std::isfinite(determinant));
                if (determinant < minimum) {
                    minimum = determinant;
                    minimum_element = element;
                }
            };
            const auto& quadrature = mfem::IntRules.Get(mfem::Geometry::CUBE, 2 * transformation->Order() + 2);
            for (int q = 0; q < quadrature.GetNPoints(); ++q) inspect(quadrature.IntPoint(q));
            for (int i = 0; i <= 2; ++i) {
                for (int j = 0; j <= 2; ++j) {
                    for (int k = 0; k <= 2; ++k) {
                        mfem::IntegrationPoint point;
                        point.Set3(i / 2.0, j / 2.0, k / 2.0);
                        inspect(point);
                    }
                }
            }
        }
        EXPECT_GT(minimum, 0.0) << "Element " << minimum_element;
    }

    double StellarVolume(stroid::StroidMesh& generated) {
        double volume = 0.0;
        for (int element = 0; element < generated.mesh->GetNE(); ++element) {
            if (generated.mesh->GetAttribute(element) == static_cast<int>(generated.config.vacuum_id.value())) continue;
            auto* transformation = generated.mesh->GetElementTransformation(element);
            const auto& quadrature = mfem::IntRules.Get(mfem::Geometry::CUBE, 3 * transformation->Order() + 3);
            for (int q = 0; q < quadrature.GetNPoints(); ++q) {
                const auto& point = quadrature.IntPoint(q);
                transformation->SetIntPoint(&point);
                volume += point.weight * transformation->Jacobian().Det();
            }
        }
        return volume;
    }

    double SurfaceRadiusError(stroid::StroidMesh& generated) {
        double error = 0.0;
        mfem::Vector physical(3);
        for (int boundary = 0; boundary < generated.mesh->GetNBE(); ++boundary) {
            if (generated.mesh->GetBdrAttribute(boundary) != static_cast<int>(generated.config.surface_bdr_id.value())) continue;
            auto* transformation = generated.mesh->GetBdrElementTransformation(boundary);
            const auto& quadrature = mfem::IntRules.Get(mfem::Geometry::SQUARE, 8);
            for (int q = 0; q < quadrature.GetNPoints(); ++q) {
                transformation->Transform(quadrature.IntPoint(q), physical);
                physical(2) /= 1.0 - generated.config.flattening.value();
                error = std::max(error, std::abs(physical.Norml2() - generated.config.r_star.value()));
            }
        }
        return error;
    }
}

TEST(NonconformingMesh, UnspecifiedVacuumLevelsPreserveUniformGeneration) {
    auto config = Configuration(2, 1);
    config.vacuum_refinement_levels.reset();
    config.vacuum_outer_refinement_levels.reset();
    auto generated = stroid::GenerateMesh(config);
    EXPECT_TRUE(generated.mesh->Conforming());
    EXPECT_EQ(generated.mesh->GetNE(), 19 * 8);
    ExpectGeometryAndCoordinateTraces(generated, false);

    config.vacuum_refinement_levels = 1;
    config.vacuum_outer_refinement_levels = 1;
    auto explicit_levels = stroid::GenerateMesh(config);
    EXPECT_EQ(explicit_levels.mesh->GetNE(), generated.mesh->GetNE());
    mfem::H1_FECollection collection(2, 3);
    mfem::FiniteElementSpace space(explicit_levels.mesh.get(), &collection);
    EXPECT_EQ(space.GetTrueVSize(), space.GetVSize());
    EXPECT_NEAR(StellarVolume(explicit_levels), StellarVolume(generated), 1.0e-12);
}

TEST(NonconformingMesh, RejectsNegativeRefinementTargets) {
    for (int field = 0; field < 3; ++field) {
        auto config = Configuration();
        if (field == 0) config.refinement_levels = -1;
        if (field == 1) config.vacuum_refinement_levels = -1;
        if (field == 2) config.vacuum_outer_refinement_levels = -1;
        EXPECT_THROW(stroid::GenerateMesh(config), std::invalid_argument);
    }
}

TEST(NonconformingMesh, DefaultOuterLevelProtectsBothSurfacesAndSavesBulkDofs) {
    auto config = Configuration(2, 3);
    config.vacuum_outer_refinement_levels.reset();
    auto generated = stroid::GenerateMesh(config);
    ASSERT_NE(generated.reference_mesh->ncmesh, nullptr);
    const auto attributes = CountAttributes(*generated.mesh);
    EXPECT_EQ(attributes.at(1), 7 * 512);
    EXPECT_EQ(attributes.at(2), 6 * 512);
    EXPECT_LT(attributes.at(3), 6 * 512);
    const auto boundaries = CountAttributes(*generated.mesh, true);
    EXPECT_EQ(boundaries.at(1), 6 * 64);
    EXPECT_EQ(boundaries.at(2), 6 * 64);
    int vacuum_minimum = 3;
    int vacuum_maximum = 0;
    for (int element = 0; element < generated.mesh->GetNE(); ++element) {
        ASSERT_EQ(generated.mesh->GetAttribute(element), generated.reference_mesh->GetAttribute(element));
        if (generated.mesh->GetAttribute(element) == 3) {
            const int depth = generated.reference_mesh->ncmesh->GetElementDepth(element);
            vacuum_minimum = std::min(vacuum_minimum, depth);
            vacuum_maximum = std::max(vacuum_maximum, depth);
        }
    }
    EXPECT_LT(vacuum_minimum, 3);
    EXPECT_EQ(vacuum_maximum, 3);
    for (int face = 0; face < generated.reference_mesh->GetNumFaces(); ++face) {
        const auto information = generated.reference_mesh->GetFaceInformation(face);
        if (!information.IsLocal()) continue;
        const int first = generated.reference_mesh->ncmesh->GetElementDepth(information.element[0].index);
        const int second = generated.reference_mesh->ncmesh->GetElementDepth(information.element[1].index);
        EXPECT_LE(std::abs(first - second), 1);
    }
    mfem::H1_FECollection collection(2, 3);
    mfem::FiniteElementSpace reduced(generated.mesh.get(), &collection);
    EXPECT_LT(reduced.GetTrueVSize(), reduced.GetVSize());
    const auto stats = stroid::stats::ComputeMeshStats(generated,
        stroid::stats::MeshStatFeatures::REFINEMENT | stroid::stats::MeshStatFeatures::CONFORMITY |
        stroid::stats::MeshStatFeatures::ELEMENT_COUNT);
    ASSERT_TRUE(stats.refinement.has_value());
    ASSERT_TRUE(stats.conformity.has_value());
    ASSERT_TRUE(stats.element_counts.has_value());
    EXPECT_EQ(stats.refinement->vacuum.min_depth, vacuum_minimum);
    EXPECT_EQ(stats.refinement->vacuum.max_depth, vacuum_maximum);
    EXPECT_EQ(stats.refinement->core.min_depth, 3);
    EXPECT_EQ(stats.refinement->core.max_depth, 3);
    EXPECT_EQ(stats.refinement->geometry_dofs, reduced.GetVSize());
    EXPECT_EQ(stats.refinement->geometry_true_dofs, reduced.GetTrueVSize());
    EXPECT_TRUE(stats.conformity->hierarchy_enabled);
    EXPECT_FALSE(stats.conformity->conforming);
    EXPECT_GT(stats.conformity->n_nonconforming_faces, 0);
    EXPECT_EQ(stats.element_counts->vacuum, attributes.at(3));
    config.vacuum_refinement_levels.reset();
    auto uniform = stroid::GenerateMesh(config);
    mfem::FiniteElementSpace full(uniform.mesh.get(), &collection);
    EXPECT_LT(reduced.GetTrueVSize(), full.GetTrueVSize());
    EXPECT_NEAR(StellarVolume(generated), StellarVolume(uniform), 2.0e-12);
    EXPECT_NEAR(SurfaceRadiusError(generated), SurfaceRadiusError(uniform), 2.0e-13);
    ExpectGeometryAndCoordinateTraces(generated);
    ExpectPositiveJacobians(*generated.mesh);
}

TEST(NonconformingMesh, OuterTargetCanExceedStellarTarget) {
    auto config = Configuration(2, 1, 0, 3);
    auto generated = stroid::GenerateMesh(config);
    EXPECT_EQ(CountAttributes(*generated.mesh, true).at(2), 6 * 64);
    ExpectGeometryAndCoordinateTraces(generated);
    ExpectPositiveJacobians(*generated.mesh);
}

TEST(NonconformingMesh, InterfaceClosureRaisesCoarseStellarBoundaryToMatchVacuum) {
    auto generated = stroid::GenerateMesh(Configuration(2, 0, 0, 3));
    ASSERT_NE(generated.reference_mesh->ncmesh, nullptr);
    EXPECT_EQ(CountAttributes(*generated.mesh, true).at(2), 6 * 64);
    int envelope_maximum = 0;
    for (int element = 0; element < generated.reference_mesh->GetNE(); ++element) {
        if (generated.reference_mesh->GetAttribute(element) == static_cast<int>(generated.config.envelope_id.value())) {
            envelope_maximum = std::max(envelope_maximum, generated.reference_mesh->ncmesh->GetElementDepth(element));
        }
    }
    EXPECT_GT(envelope_maximum, 0);
    for (int face = 0; face < generated.reference_mesh->GetNumFaces(); ++face) {
        const auto information = generated.reference_mesh->GetFaceInformation(face);
        if (!information.IsLocal()) continue;
        const int first = generated.reference_mesh->ncmesh->GetElementDepth(information.element[0].index);
        const int second = generated.reference_mesh->ncmesh->GetElementDepth(information.element[1].index);
        EXPECT_LE(std::abs(first - second), 1);
    }
    ExpectGeometryAndCoordinateTraces(generated);
    ExpectPositiveJacobians(*generated.mesh);
}

TEST(NonconformingMesh, RefinementAndExteriorCoordinateAreInvariantUnderSmallLengthScales) {
    auto config = Configuration(2, 1, 0, 3);
    auto reference = stroid::GenerateMesh(config);
    constexpr double scale = 1.0e-15;
    config.r_core = config.r_core.value() * scale;
    config.r_star = config.r_star.value() * scale;
    config.r_infinity = config.r_infinity.value() * scale;
    auto scaled = stroid::GenerateMesh(config);
    ASSERT_EQ(scaled.mesh->GetNE(), reference.mesh->GetNE());
    ASSERT_EQ(scaled.mesh->GetNodes()->Size(), reference.mesh->GetNodes()->Size());
    double coordinate_error = 0.0;
    for (int dof = 0; dof < scaled.mesh->GetNodes()->Size(); ++dof) {
        coordinate_error = std::max(coordinate_error,
            std::abs((*scaled.mesh->GetNodes())(dof) / scale - (*reference.mesh->GetNodes())(dof)));
    }
    EXPECT_LT(coordinate_error, 5.0e-12);
    ASSERT_EQ(scaled.exterior_coordinate->values->Size(), reference.exterior_coordinate->values->Size());
    mfem::Vector difference(*scaled.exterior_coordinate->values);
    difference -= *reference.exterior_coordinate->values;
    EXPECT_LT(difference.Normlinf(), 5.0e-12);
    ExpectGeometryAndCoordinateTraces(scaled);
    ExpectPositiveJacobians(*scaled.mesh);
}

TEST(NonconformingMesh, CurvedGeometryAndScalarConstraintsAcrossOrdersAndMappings) {
    for (const std::string mapping : {"multi_block", "spherified"}) {
        for (const int order : {1, 2, 3, 4}) {
            SCOPED_TRACE(mapping + " order=" + std::to_string(order));
            auto config = Configuration(order);
            config.core_mapping = mapping;
            config.flattening = 0.2;
            config.core_id = 11;
            config.envelope_id = 17;
            config.vacuum_id = 23;
            config.surface_bdr_id = 31;
            config.inf_bdr_id = 37;
            auto generated = stroid::GenerateMesh(config);
            EXPECT_EQ(CountAttributes(*generated.mesh).size(), 3);
            EXPECT_EQ(CountAttributes(*generated.mesh, true).size(), 2);
            ExpectGeometryAndCoordinateTraces(generated);
            // The legacy spherified core has known corner degeneracies. Its
            // envelope and vacuum must still remain strictly oriented.
            ExpectPositiveJacobians(*generated.mesh, mapping == "spherified" ? 11 : -1);
        }
    }
}

TEST(NonconformingMesh, NoVacuumLeavesOnlyTheUniformStellarMesh) {
    auto config = Configuration(2, 1, 0, 3);
    config.include_external_domain = false;
    EXPECT_THROW(stroid::GenerateMesh(config), std::invalid_argument);
    config.vacuum_refinement_levels.reset();
    config.vacuum_outer_refinement_levels.reset();
    auto generated = stroid::GenerateMesh(config);
    EXPECT_EQ(generated.mesh->GetNE(), 13 * 8);
    EXPECT_EQ(generated.exterior_coordinate, nullptr);
    EXPECT_EQ(CountAttributes(*generated.mesh).size(), 2);
    ExpectPositiveJacobians(*generated.mesh);
}

TEST(NonconformingMesh, SerializationAndSubsequentRefinementPreserveHierarchy) {
    auto generated = stroid::GenerateMesh(Configuration());
    const auto path = std::filesystem::temp_directory_path() / "stroid_nonconforming_round_trip.smesh";
    stroid::IO::SaveStroidMesh(generated, path.string(), "Nonconforming hierarchy regression");
    auto result = stroid::IO::LoadStroidMesh(path.string());
    ASSERT_TRUE(result.has_value()) << result.error();
    auto loaded = std::move(*result);
    ASSERT_NE(loaded.reference_mesh->ncmesh, nullptr);
    ASSERT_NE(loaded.mesh->ncmesh, nullptr);
    EXPECT_EQ(loaded.config.vacuum_refinement_levels, generated.config.vacuum_refinement_levels);
    EXPECT_EQ(loaded.config.vacuum_outer_refinement_levels, generated.config.vacuum_outer_refinement_levels);
    ASSERT_EQ(loaded.mesh->GetNE(), generated.mesh->GetNE());
    for (int element = 0; element < loaded.mesh->GetNE(); ++element) {
        EXPECT_EQ(loaded.reference_mesh->ncmesh->GetElementDepth(element),
                  generated.reference_mesh->ncmesh->GetElementDepth(element));
    }
    EXPECT_NEAR(StellarVolume(loaded), StellarVolume(generated), 2.0e-12);
    ExpectGeometryAndCoordinateTraces(loaded);
    stroid::refinement::UniformRefinement(loaded, 1);
    EXPECT_EQ(loaded.mesh->GetNE(), generated.mesh->GetNE() * 8);
    EXPECT_EQ(loaded.mesh->GetNE(), loaded.reference_mesh->GetNE());
    ExpectGeometryAndCoordinateTraces(loaded);
    ExpectPositiveJacobians(*loaded.mesh);
    std::error_code error;
    std::filesystem::remove(path, error);
    EXPECT_FALSE(error);
}

TEST(NonconformingMesh, LinearPhysicalPatchSolveUsesIndependentDofs) {
    auto generated = stroid::GenerateMesh(Configuration());
    mfem::H1_FECollection collection(2, 3);
    mfem::FiniteElementSpace space(generated.mesh.get(), &collection);
    ASSERT_LT(space.GetTrueVSize(), space.GetVSize());
    mfem::FunctionCoefficient exact([](const mfem::Vector& point) {
        return 1.0 + 0.3 * point(0) - 0.2 * point(1) + 0.1 * point(2);
    });
    mfem::Array<int> boundary(generated.mesh->bdr_attributes.Max());
    boundary = 0;
    boundary[static_cast<int>(generated.config.inf_bdr_id.value()) - 1] = 1;
    mfem::Array<int> essential;
    space.GetEssentialTrueDofs(boundary, essential);
    mfem::GridFunction solution(&space);
    solution = 0.0;
    solution.ProjectBdrCoefficient(exact, boundary);
    mfem::LinearForm rhs(&space);
    rhs = 0.0;
    mfem::ConstantCoefficient one(1.0);
    mfem::BilinearForm form(&space);
    const auto& quadrature = mfem::IntRules.Get(mfem::Geometry::CUBE, 10);
    auto* diffusion = new mfem::DiffusionIntegrator(one);
    diffusion->SetIntRule(&quadrature);
    form.AddDomainIntegrator(diffusion);
    form.Assemble();
    mfem::OperatorPtr system;
    mfem::Vector independent, forcing;
    form.FormLinearSystem(essential, solution, rhs, system, independent, forcing);
    EXPECT_EQ(system->Height(), space.GetTrueVSize());
    mfem::GSSmoother preconditioner(static_cast<mfem::SparseMatrix&>(*system));
    mfem::CGSolver solver;
    solver.SetOperator(*system);
    solver.SetPreconditioner(preconditioner);
    solver.SetRelTol(1.0e-13);
    solver.SetAbsTol(1.0e-14);
    solver.SetMaxIter(1500);
    solver.SetPrintLevel(-1);
    solver.Mult(forcing, independent);
    ASSERT_TRUE(solver.GetConverged());
    mfem::Vector residual(forcing.Size());
    system->Mult(independent, residual);
    residual -= forcing;
    EXPECT_LT(residual.Norml2() / forcing.Norml2(), 2.0e-12);
    form.RecoverFEMSolution(independent, rhs, solution);
    EXPECT_LT(solution.ComputeL2Error(exact), 2.0e-8);
    ExpectConstrainedField(solution);
}

TEST(NonconformingMesh, StellarVolumeAndSurfaceShapeConverge) {
    auto coarse_config = Configuration(2, 2);
    coarse_config.vacuum_outer_refinement_levels.reset();
    auto coarse = stroid::GenerateMesh(coarse_config);
    auto fine_config = coarse_config;
    fine_config.refinement_levels = 3;
    auto fine = stroid::GenerateMesh(fine_config);
    const double exact_volume = 4.0 * kPi / 3.0;
    const double coarse_volume_error = std::abs(StellarVolume(coarse) - exact_volume);
    const double fine_volume_error = std::abs(StellarVolume(fine) - exact_volume);
    EXPECT_GT(coarse_volume_error, 1.0e-10);
    EXPECT_LT(fine_volume_error, 0.5 * coarse_volume_error);
    EXPECT_LT(SurfaceRadiusError(fine), 0.5 * SurfaceRadiusError(coarse));
}

TEST(NonconformingMesh, TMOPPreservesHangingConstraintsAndBoundaryTraces) {
    auto config = Configuration(1);
    auto initial = stroid::GenerateMesh(config);
    config.optimization_methods = stroid::config::OptimizationMethods{true, true};
    auto generated = stroid::GenerateMesh(config);
    mfem::Array<int> marker(generated.mesh->bdr_attributes.Max());
    marker = 1;
    mfem::Array<int> essential;
    generated.mesh->GetNodalFESpace()->GetEssentialTrueDofs(marker, essential);
    mfem::Vector initial_nodes, optimized_nodes;
    initial.mesh->GetNodes()->GetTrueDofs(initial_nodes);
    generated.mesh->GetNodes()->GetTrueDofs(optimized_nodes);
    ASSERT_EQ(initial_nodes.Size(), optimized_nodes.Size());
    for (const int dof : essential) EXPECT_NEAR(initial_nodes(dof), optimized_nodes(dof), 2.0e-13);
    ExpectGeometryAndCoordinateTraces(generated);
    ExpectPositiveJacobians(*generated.mesh);
}
