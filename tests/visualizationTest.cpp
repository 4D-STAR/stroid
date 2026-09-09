#include <gtest/gtest.h>

#include "stroid/stroid.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

namespace {
    std::string Serialize(const mfem::Mesh& mesh) {
        std::ostringstream stream;
        stream.precision(std::numeric_limits<double>::max_digits10);
        mesh.Print(stream);
        return stream.str();
    }

    int HangingFaces(const mfem::Mesh& mesh) {
        int count = 0;
        for (int face = 0; face < mesh.GetNumFaces(); ++face) {
            count += mesh.GetFaceInformation(face).IsNonconformingFine();
        }
        return count;
    }

    void TessellatedPoint(mfem::ElementTransformation& transformation,
                          const mfem::IntegrationPoint& point, int subdivisions,
                          mfem::Vector& value) {
        const double coordinates[] = {point.x, point.y, point.z};
        int cell[3];
        double local[3];
        for (int d = 0; d < 3; ++d) {
            const double scaled = coordinates[d] * subdivisions;
            cell[d] = std::clamp(static_cast<int>(std::floor(scaled)), 0, subdivisions - 1);
            local[d] = scaled - cell[d];
        }
        value = 0.0;
        mfem::Vector corner(3);
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                for (int k = 0; k < 2; ++k) {
                    mfem::IntegrationPoint sample;
                    sample.Set3(static_cast<double>(cell[0] + i) / subdivisions,
                                static_cast<double>(cell[1] + j) / subdivisions,
                                static_cast<double>(cell[2] + k) / subdivisions);
                    transformation.Transform(sample, corner);
                    const double weight = (i ? local[0] : 1.0 - local[0]) *
                                          (j ? local[1] : 1.0 - local[1]) *
                                          (k ? local[2] : 1.0 - local[2]);
                    value.Add(weight, corner);
                }
            }
        }
    }

    double FaceGap(mfem::Mesh& mesh, int subdivisions = 0) {
        double maximum = 0.0;
        mfem::Vector first(3), second(3);
        const auto& quadrature = mfem::IntRules.Get(mfem::Geometry::SQUARE, 6);
        for (int face = 0; face < mesh.GetNumFaces(); ++face) {
            if (!mesh.GetFaceInformation(face).IsLocal()) continue;
            auto* transformation = mesh.GetFaceElementTransformations(face);
            for (int q = 0; q < quadrature.GetNPoints(); ++q) {
                transformation->SetAllIntPoints(&quadrature.IntPoint(q));
                const auto first_point = transformation->Elem1->GetIntPoint();
                const auto second_point = transformation->Elem2->GetIntPoint();
                if (subdivisions > 0) {
                    TessellatedPoint(*transformation->Elem1, first_point, subdivisions, first);
                    TessellatedPoint(*transformation->Elem2, second_point, subdivisions, second);
                } else {
                    transformation->Elem1->Transform(first_point, first);
                    transformation->Elem2->Transform(second_point, second);
                }
                first -= second;
                maximum = std::max(maximum, first.Norml2());
            }
        }
        return maximum;
    }

    stroid::StroidMesh CurvedVacuumMesh() {
        stroid::config::MeshConfig config;
        config.order = 3;
        config.refinement_levels = 1;
        config.vacuum_refinement_levels = 0;
        config.vacuum_outer_refinement_levels = 2;
        config.flattening = 0.15;
        config.optimization_methods = stroid::config::OptimizationMethods{false, true};
        return stroid::GenerateMesh(config);
    }
}

TEST(Visualization, ConformingDisplayEliminatesCurvedTessellationGaps) {
    auto generated = CurvedVacuumMesh();
    auto& source = *generated.mesh;
    ASSERT_GT(HangingFaces(source), 0);
    EXPECT_LT(FaceGap(source), 5.0e-12);
    // Reproduce visible cracks even though the finite-element traces coincide.
    EXPECT_GT(FaceGap(source, 2), 1.0e-4);
    const auto original = Serialize(source);
    const auto original_reference = Serialize(*generated.reference_mesh);

    auto display = stroid::IO::MakeConformingVisualizationMesh(source);
    ASSERT_NE(display, nullptr);
    EXPECT_GT(display->GetNE(), source.GetNE());
    EXPECT_EQ(HangingFaces(*display), 0);
    EXPECT_EQ(display->GetNodalFESpace()->GetNDofs(),
              display->GetNodalFESpace()->GetTrueVSize() / display->SpaceDimension());
    EXPECT_LT(FaceGap(*display), 5.0e-12);
    for (int subdivisions : {1, 2, 3, 4}) {
        EXPECT_LT(FaceGap(*display, subdivisions), 5.0e-12) << subdivisions;
    }

    std::istringstream stream(Serialize(*display));
    mfem::Mesh reloaded(stream, 1, 1, true);
    EXPECT_EQ(HangingFaces(reloaded), 0);
    EXPECT_LT(FaceGap(reloaded, 2), 5.0e-12);
    EXPECT_EQ(Serialize(source), original);
    EXPECT_EQ(Serialize(*generated.reference_mesh), original_reference);
}

TEST(Visualization, DisplayRefinementRestrictsExistingGeometryAndAttributes) {
    auto generated = CurvedVacuumMesh();
    auto display = stroid::IO::MakeConformingVisualizationMesh(*generated.mesh);
    auto reference = stroid::IO::MakeConformingVisualizationMesh(*generated.reference_mesh);
    ASSERT_EQ(display->GetNE(), reference->GetNE());

    mfem::DenseMatrix centers(3, reference->GetNE());
    mfem::Vector point(3), expected(3), actual(3);
    for (int element = 0; element < reference->GetNE(); ++element) {
        reference->GetElementCenter(element, point);
        centers.SetCol(element, point);
    }
    mfem::Array<int> parents;
    mfem::Array<mfem::IntegrationPoint> parent_points;
    ASSERT_EQ(generated.reference_mesh->FindPoints(centers, parents, parent_points, false),
              reference->GetNE());
    double maximum_error = 0.0;
    for (int element = 0; element < display->GetNE(); ++element) {
        const int parent = parents[element];
        ASSERT_GE(parent, 0);
        EXPECT_EQ(display->GetAttribute(element), generated.mesh->GetAttribute(parent));
        mfem::InverseElementTransformation inverse(
            generated.reference_mesh->GetElementTransformation(parent));
        auto* logical = reference->GetElementTransformation(element);
        auto* physical = display->GetElementTransformation(element);
        auto* original = generated.mesh->GetElementTransformation(parent);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    mfem::IntegrationPoint sample, parent_sample;
                    sample.Set3(i / 2.0, j / 2.0, k / 2.0);
                    logical->Transform(sample, point);
                    ASSERT_EQ(inverse.Transform(point, parent_sample),
                              mfem::InverseElementTransformation::Inside);
                    original->Transform(parent_sample, expected);
                    physical->Transform(sample, actual);
                    actual -= expected;
                    maximum_error = std::max(maximum_error, actual.Norml2());
                }
            }
        }
    }
    EXPECT_LT(maximum_error, 5.0e-12);
}

TEST(Visualization, AlreadyMatchingFacesNeedNoAdditionalElements) {
    stroid::config::MeshConfig config;
    config.order = 3;
    config.refinement_levels = 1;
    config.optimization_methods = stroid::config::OptimizationMethods{false, true};
    for (bool hierarchy : {false, true}) {
        config.vacuum_refinement_levels = hierarchy ? std::optional<int>(1) : std::nullopt;
        auto generated = stroid::GenerateMesh(config);
        if (hierarchy) generated.mesh->EnsureNCMesh();
        const auto original = Serialize(*generated.mesh);
        auto display = stroid::IO::MakeConformingVisualizationMesh(*generated.mesh);
        EXPECT_EQ(display->GetNE(), generated.mesh->GetNE());
        EXPECT_EQ(Serialize(*display), original);
        EXPECT_EQ(Serialize(*generated.mesh), original);
    }
}
