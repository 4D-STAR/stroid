#include <gtest/gtest.h>

#include "fourdst/config/config.h"
#include "stroid/config/config.h"
#include "stroid/IO/mesh.h"
#include "stroid/topology/curvilinear.h"
#include "stroid/topology/mapping.h"
#include "stroid/topology/topology.h"
#include "stroid/utils/mesh_utils.h"
#include "stroid/stroid.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <array>
#include <limits>
#include <print>
#include <complex>
#include <stroid/stroid.h>

#include "stroid/utils/mesh_stats.h"
#include "stroid/utils/types.h"

namespace {

constexpr double kPi = 3.14159265358979323846;

using Config = fourdst::config::Config<stroid::config::MeshConfig>;


std::filesystem::path GetSourceRoot() {
    if (const char* env = std::getenv("MESON_SOURCE_ROOT")) {
        return {env};
    }
    return std::filesystem::current_path();
}

std::unique_ptr<Config> LoadConfigFromRepo(const std::filesystem::path& relative_path) {
    auto cfg_ptr = std::make_unique<Config>();
    cfg_ptr->load((GetSourceRoot() / relative_path).string());
    return cfg_ptr;
}


bool IsFiniteMeshNodes(const mfem::Mesh& mesh) {
    const mfem::GridFunction* nodes = mesh.GetNodes();
    if (!nodes) {
        return false;
    }
    const int vdim = nodes->FESpace()->GetVDim();
    const int ndofs = nodes->FESpace()->GetNDofs();
    for (int i = 0; i < ndofs; ++i) {
        for (int d = 0; d < vdim; ++d) {
            const double val = (*nodes)(nodes->FESpace()->DofToVDof(i, d));
            if (!std::isfinite(val)) {
                return false;
            }
        }
    }
    return true;
}

std::map<int, int> CountVolumeAttributes(const mfem::Mesh& mesh) {
    std::map<int, int> counts;
    for (int i = 0; i < mesh.GetNE(); ++i) {
        counts[mesh.GetAttribute(i)]++;
    }
    return counts;
}

std::map<int, int> CountBoundaryAttributes(const mfem::Mesh& mesh) {
    std::map<int, int> counts;
    for (int i = 0; i < mesh.GetNBE(); ++i) {
        counts[mesh.GetBdrAttribute(i)]++;
    }
    return counts;
}

mfem::Vector TransformCopy(const mfem::Vector& in, const Config& cfg, int attribute_id = 0) {
    mfem::Vector out = in;
    stroid::topology::TransformPoint(out, cfg, attribute_id);
    return out;
}

double IntegrateElementVolume(const mfem::Mesh& mesh, int element_id) {
    mfem::ElementTransformation* T = const_cast<mfem::Mesh&>(mesh).GetElementTransformation(element_id);
    const int order = std::max(2, 2 * T->Order() + 2);
    const mfem::IntegrationRule& ir = mfem::IntRules.Get(T->GetGeometryType(), order);

    double volume = 0.0;
    for (int j = 0; j < ir.GetNPoints(); ++j) {
        const mfem::IntegrationPoint& ip = ir.IntPoint(j);
        T->SetIntPoint(&ip);
        volume += ip.weight * std::abs(T->Weight());
    }
    return volume;
}

double ComputeMeshVolume(const mfem::Mesh& mesh) {
    double volume = 0.0;
    for (int i = 0; i < mesh.GetNE(); ++i) {
        volume += IntegrateElementVolume(mesh, i);
    }
    return volume;
}

double ComputeMeshVolumeForAttributes(const mfem::Mesh& mesh, const std::set<int>& attributes) {
    double volume = 0.0;
    for (int i = 0; i < mesh.GetNE(); ++i) {
        if (!attributes.contains(mesh.GetAttribute(i))) {
            continue;
        }
        volume += IntegrateElementVolume(mesh, i);
    }
    return volume;
}

std::unique_ptr<mfem::Mesh> BuildProjectedMesh(const Config& cfg) {
    std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    stroid::topology::Finalize(*mesh, cfg);
    stroid::topology::PromoteToHighOrder(*mesh, cfg);
    stroid::topology::ProjectMesh(*mesh, cfg);
    return mesh;
}

void ExpectExteriorCoordinateRange(stroid::StroidMesh& stroid_mesh) {
    ASSERT_NE(stroid_mesh.mesh, nullptr);
    ASSERT_NE(stroid_mesh.exterior_coordinate, nullptr);
    ASSERT_NE(stroid_mesh.exterior_coordinate->space, nullptr);
    ASSERT_NE(stroid_mesh.exterior_coordinate->values, nullptr);
    ASSERT_EQ(stroid_mesh.exterior_coordinate->space->GetMesh(), stroid_mesh.mesh.get());
    ASSERT_EQ(stroid_mesh.exterior_coordinate->values->FESpace(), stroid_mesh.exterior_coordinate->space.get());

    mfem::Mesh& mesh = *stroid_mesh.mesh;
    mfem::GridFunction& coordinate = *stroid_mesh.exterior_coordinate->values;
    const int vacuum_attribute = static_cast<int>(stroid_mesh.config.vacuum_id.value());
    bool sampled_vacuum = false;

    for (int element_id = 0; element_id < mesh.GetNE(); ++element_id) {
        const mfem::FiniteElement& element = *stroid_mesh.exterior_coordinate->space->GetFE(element_id);
        const mfem::IntegrationRule& integration_rule = mfem::IntRules.Get(element.GetGeomType(), 2 * element.GetOrder() + 4);

        for (int q = 0; q < integration_rule.GetNPoints(); ++q) {
            const double value = coordinate.GetValue(element_id, integration_rule.IntPoint(q));
            EXPECT_TRUE(std::isfinite(value));

            if (mesh.GetAttribute(element_id) == vacuum_attribute) {
                sampled_vacuum = true;
                EXPECT_GE(value, -1.0e-12);
                EXPECT_LE(value, 1.0 + 1.0e-12);
            } else {
                EXPECT_NEAR(value, 0.0, 1.0e-12);
            }
        }
    }

    EXPECT_TRUE(sampled_vacuum);
}

void ExpectExteriorCoordinateBoundaryTraces(stroid::StroidMesh& stroid_mesh) {
    ASSERT_NE(stroid_mesh.mesh, nullptr);
    ASSERT_NE(stroid_mesh.exterior_coordinate, nullptr);
    ASSERT_NE(stroid_mesh.exterior_coordinate->values, nullptr);

    mfem::Mesh& mesh = *stroid_mesh.mesh;
    mfem::GridFunction& coordinate = *stroid_mesh.exterior_coordinate->values;
    const int vacuum_attribute = static_cast<int>(stroid_mesh.config.vacuum_id.value());
    const int infinity_boundary = static_cast<int>(stroid_mesh.config.inf_bdr_id.value());
    int stellar_vacuum_faces = 0;
    int infinity_faces = 0;

    for (int face_id = 0; face_id < mesh.GetNumFaces(); ++face_id) {
        mfem::FaceElementTransformations* transformation = mesh.GetFaceElementTransformations(face_id);
        if (transformation == nullptr || transformation->Elem1 == nullptr || transformation->Elem2 == nullptr) continue;

        const bool element_1_vacuum = transformation->Elem1->Attribute == vacuum_attribute;
        const bool element_2_vacuum = transformation->Elem2->Attribute == vacuum_attribute;
        if (element_1_vacuum == element_2_vacuum) continue;

        ++stellar_vacuum_faces;
        const mfem::IntegrationRule& integration_rule = mfem::IntRules.Get(transformation->GetGeometryType(), 6);

        for (int q = 0; q < integration_rule.GetNPoints(); ++q) {
            const mfem::IntegrationPoint& face_point = integration_rule.IntPoint(q);
            transformation->SetAllIntPoints(&face_point);

            const int vacuum_element = element_1_vacuum ? transformation->Elem1No : transformation->Elem2No;
            const mfem::IntegrationPoint& vacuum_point = element_1_vacuum ? transformation->Elem1->GetIntPoint() : transformation->Elem2->GetIntPoint();
            EXPECT_NEAR(coordinate.GetValue(vacuum_element, vacuum_point), 0.0, 1.0e-12);
        }
    }

    for (int boundary_element = 0; boundary_element < mesh.GetNBE(); ++boundary_element) {
        if (mesh.GetBdrAttribute(boundary_element) != infinity_boundary) continue;

        mfem::FaceElementTransformations* transformation = mesh.GetBdrFaceTransformations(boundary_element);
        ASSERT_NE(transformation, nullptr);
        ASSERT_NE(transformation->Elem1, nullptr);
        ++infinity_faces;

        const mfem::IntegrationRule& integration_rule = mfem::IntRules.Get(transformation->GetGeometryType(), 6);

        for (int q = 0; q < integration_rule.GetNPoints(); ++q) {
            const mfem::IntegrationPoint& face_point = integration_rule.IntPoint(q);
            transformation->SetAllIntPoints(&face_point);
            EXPECT_NEAR(coordinate.GetValue(transformation->Elem1No, transformation->Elem1->GetIntPoint()), 1.0, 1.0e-12);
        }
    }

    EXPECT_GT(stellar_vacuum_faces, 0);
    EXPECT_GT(infinity_faces, 0);
}

double ComputeStellarVolumeWithDomainLFIntegrator(mfem::Mesh& mesh, const Config& cfg) {
    const int mesh_max_attr = mesh.attributes.Size() > 0 ? mesh.attributes.Max() : 0;
    const int cfg_max_attr = static_cast<int>(std::max({cfg->core_id.value(), cfg->envelope_id.value(), cfg->vacuum_id.value()}));
    const int coeff_size = std::max(1, std::max(mesh_max_attr, cfg_max_attr));

    mfem::Vector attr_coeff(coeff_size);
    attr_coeff = 0.0;
    attr_coeff(static_cast<int>(cfg->core_id.value()) - 1) = 1.0;
    attr_coeff(static_cast<int>(cfg->envelope_id.value()) - 1) = 1.0;

    mfem::PWConstCoefficient stellar_coeff(attr_coeff);
    mfem::L2_FECollection fec(0, mesh.Dimension());
    mfem::FiniteElementSpace fes(&mesh, &fec);
    mfem::LinearForm lf(&fes);
    lf.AddDomainIntegrator(new mfem::DomainLFIntegrator(stellar_coeff));
    lf.Assemble();
    return lf.Sum();
}

double ColumnNorm2(const mfem::DenseMatrix& J, int col) {
    double n2 = 0.0;
    for (int r = 0; r < J.Height(); ++r) {
        n2 += J(r, col) * J(r, col);
    }
    return std::sqrt(n2);
}

double ComputeHexEdgeRatio(const mfem::Mesh& mesh, int elem_id) {
    static constexpr std::array<std::array<int, 2>, 12> edges = {{
        {{0, 1}}, {{1, 2}}, {{2, 3}}, {{3, 0}},
        {{4, 5}}, {{5, 6}}, {{6, 7}}, {{7, 4}},
        {{0, 4}}, {{1, 5}}, {{2, 6}}, {{3, 7}}
    }};

    const mfem::Element* e = mesh.GetElement(elem_id);
    if (e->GetNVertices() != 8) {
        return 1.0;
    }

    const int* v = e->GetVertices();
    double min_edge = std::numeric_limits<double>::infinity();
    double max_edge = 0.0;

    for (const auto& [a, b] : edges) {
        const double* pa = mesh.GetVertex(v[a]);
        const double* pb = mesh.GetVertex(v[b]);
        const double dx = pa[0] - pb[0];
        const double dy = pa[1] - pb[1];
        const double dz = pa[2] - pb[2];
        const double len = std::sqrt(dx * dx + dy * dy + dz * dz);
        min_edge = std::min(min_edge, len);
        max_edge = std::max(max_edge, len);
    }

    if (min_edge <= 0.0 || !std::isfinite(min_edge)) {
        return std::numeric_limits<double>::infinity();
    }
    return max_edge / min_edge;
}

struct ConditioningStats {
    double min_det = std::numeric_limits<double>::infinity();
    double max_det = 0.0;
    double min_scaled_jac = std::numeric_limits<double>::infinity();
    double max_stretch_ratio = 0.0;
    double max_edge_ratio = 0.0;
    int samples = 0;
};

ConditioningStats CollectConditioningStats(const mfem::Mesh& mesh, const std::set<int>& attrs) {
    ConditioningStats stats;

    for (int i = 0; i < mesh.GetNE(); ++i) {
        if (!attrs.empty() && !attrs.contains(mesh.GetAttribute(i))) {
            continue;
        }

        stats.max_edge_ratio = std::max(stats.max_edge_ratio, ComputeHexEdgeRatio(mesh, i));

        mfem::ElementTransformation* T = const_cast<mfem::Mesh&>(mesh).GetElementTransformation(i);
        const int order = std::max(2, 2 * T->Order() + 2);
        const mfem::IntegrationRule& ir = mfem::IntRules.Get(T->GetGeometryType(), order);

        for (int j = 0; j < ir.GetNPoints(); ++j) {
            const mfem::IntegrationPoint& ip = ir.IntPoint(j);
            T->SetIntPoint(&ip);

            const mfem::DenseMatrix& J = T->Jacobian();
            const double det = T->Weight();
            const double abs_det = std::abs(det);

            const double c0 = ColumnNorm2(J, 0);
            const double c1 = ColumnNorm2(J, 1);
            const double c2 = ColumnNorm2(J, 2);

            const double denom = c0 * c1 * c2;
            const double scaled_jac = (denom > 0.0) ? (abs_det / denom) : 0.0;

            const double cmax = std::max({c0, c1, c2});
            const double cmin = std::max(1e-16, std::min({c0, c1, c2}));
            const double stretch_ratio = cmax / cmin;

            stats.min_det = std::min(stats.min_det, det);
            stats.max_det = std::max(stats.max_det, abs_det);
            stats.min_scaled_jac = std::min(stats.min_scaled_jac, scaled_jac);
            stats.max_stretch_ratio = std::max(stats.max_stretch_ratio, stretch_ratio);
            stats.samples++;
        }
    }

    return stats;
}

std::optional<double> EvalGridFunctionAtPoint(
    mfem::Mesh& mesh,
    const mfem::Vector& x,
    const mfem::GridFunction& u ){

    mfem::Array<int> elem_ids;
    mfem::Array<mfem::IntegrationPoint> ips;
    mfem::DenseMatrix P(x.Size(), 1);
    P.SetCol(0, x);

    mesh.FindPoints(P, elem_ids, ips, false);

    if (elem_ids.Size() > 0 && elem_ids[0] >= 0) {
        return u.GetValue(elem_ids[0], ips[0]);
    } else {
        return std::nullopt;

    }
}

} // namespace

/**
 * @brief Test suite for the Stroid library
 */
class stroidTest : public ::testing::Test {};

/**
 * @brief Verifies the default multi-block topology including the vacuum.
 * @details
 * Rationale: this is the fastest canary for accidental edits in block construction order,
 * vertex indexing, or boundary-face assembly.
 * Method: build the default skeleton and assert exact counts (3D, 32 vertices, 19 hexes, 12 bdr quads).
 * If this fails: inspect `stroid::topology::BuildSkeleton` in `src/lib/topology/topology.cpp`,
 * especially `add_box`, `stellar_shells`, and `surface_bdr_quads`, plus ID defaults in
 * `src/include/stroid/config/config.h`.
 */
TEST_F(stroidTest, BuildSkeleton_DefaultCounts) {
    const Config cfg;
    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);

    ASSERT_NE(mesh, nullptr);
    EXPECT_EQ(mesh->Dimension(), 3);
    EXPECT_EQ(mesh->GetNV(), 32);
    EXPECT_EQ(mesh->GetNE(), 19);
    EXPECT_EQ(mesh->GetNBE(), 12);
}

/**
 * @brief Verifies topology cardinalities when the external vacuum domain is enabled.
 * @details
 * Rationale: external-domain regressions usually surface first as wrong element/boundary counts.
 * Method: load `configs/test_external_domain.toml`, build skeleton, assert exact counts (24, 13, 12).
 * If this fails: inspect vacuum block creation and infinity boundary insertion in
 * `src/lib/topology/topology.cpp` (`vacuum_shells`, `inf_bdr_quads`) and config parsing path.
 */
TEST_F(stroidTest, BuildSkeleton_ExternalDomainCounts) {
    const auto cfg_ptr = LoadConfigFromRepo("configs/test_external_domain.toml");
    const auto& cfg = *cfg_ptr;
    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);

    ASSERT_NE(mesh, nullptr);
    EXPECT_EQ(mesh->Dimension(), 3);
    EXPECT_EQ(mesh->GetNV(), 24);
    EXPECT_EQ(mesh->GetNE(), 13);
    EXPECT_EQ(mesh->GetNBE(), 12);
}

/**
 * @brief Confirms attribute bookkeeping for core/envelope/vacuum and surface/infinity boundaries.
 * @details
 * Rationale: physics coupling depends on stable material and boundary IDs, not just geometry.
 * Method: count attributes immediately after skeleton build and assert expected multiplicities.
 * If this fails: inspect element insertion attribute arguments in `BuildSkeleton` and verify
 * `core_id`, `envelope_id`, `vacuum_id`, `surface_bdr_id`, `inf_bdr_id` in config fixtures.
 */
TEST_F(stroidTest, BuildSkeleton_ExternalDomainAttributes) {
    const auto cfg_ptr = LoadConfigFromRepo("configs/test_external_domain.toml");
    const auto& cfg = *cfg_ptr;
    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);

    ASSERT_NE(mesh, nullptr);

    const auto volume_attr_counts = CountVolumeAttributes(*mesh);
    EXPECT_EQ(volume_attr_counts.at(static_cast<int>(cfg->core_id.value())), 1);
    EXPECT_EQ(volume_attr_counts.at(static_cast<int>(cfg->envelope_id.value())), 6);
    EXPECT_EQ(volume_attr_counts.at(static_cast<int>(cfg->vacuum_id.value())), 6);

    const auto boundary_attr_counts = CountBoundaryAttributes(*mesh);
    EXPECT_EQ(boundary_attr_counts.at(static_cast<int>(cfg->surface_bdr_id.value())), 6);
    EXPECT_EQ(boundary_attr_counts.at(static_cast<int>(cfg->inf_bdr_id.value())), 6);
}


/**
 * @brief Ensures `Finalize` performs refinement and preserves conforming topology.
 * @details
 * Rationale: `Finalize` is the topology gate before high-order projection; nonconforming output here
 * contaminates every downstream stage.
 * Method: compare element count pre/post finalize and assert `mesh.Conforming()`.
 * If this fails: inspect `stroid::topology::Finalize` in `src/lib/topology/topology.cpp`
 * (`FinalizeTopology`, orientation checks, `UniformRefinement` loop).
 */
TEST_F(stroidTest, Finalize_RefinementIncreasesElements) {
    const Config cfg;
    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    const int initial_elements = mesh->GetNE();

    stroid::topology::Finalize(*mesh, cfg);

    EXPECT_GT(mesh->GetNE(), initial_elements);
    EXPECT_TRUE(mesh->Conforming());
}

/**
 * @brief Checks exact 3D uniform-refinement scaling for default topology.
 * @details
 * Rationale: each hexahedron should split into 8; this catches subtle refine-loop regressions.
 * Method: run with fixed `refinement_levels=2` config and assert `NE_final = NE_initial * 8^2`.
 * If this fails: inspect refine-loop count and any topology-side early exits in `Finalize`.
 */
TEST_F(stroidTest, Finalize_DefaultRefinementScalesHexCountByEightPowerL) {
    const auto cfg_ptr = LoadConfigFromRepo("configs/test_refinement_l2.toml");
    const auto& cfg = *cfg_ptr;

    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    const int initial_elements = mesh->GetNE();

    stroid::topology::Finalize(*mesh, cfg);

    const int expected = initial_elements * 8 * 8;
    EXPECT_EQ(mesh->GetNE(), expected);
}

/**
 * @brief Checks exact 3D uniform-refinement scaling for external-domain topology.
 * @details
 * Rationale: refinement behavior must be independent of whether vacuum blocks are present.
 * Method: use fixed `refinement_levels=1` external config and assert `NE_final = NE_initial * 8`.
 * If this fails: inspect `Finalize` and verify external-domain elements are not excluded from refinement.
 */
TEST_F(stroidTest, Finalize_ExternalDomainRefinementScalesHexCountByEightPowerL) {
    const auto cfg_ptr = LoadConfigFromRepo("configs/test_external_domain_refinement_l1.toml");
    const auto& cfg = *cfg_ptr;

    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    const int initial_elements = mesh->GetNE();

    stroid::topology::Finalize(*mesh, cfg);

    const int expected = initial_elements * 8;
    EXPECT_EQ(mesh->GetNE(), expected);
}

/**
 * @brief Validates conformance + attribute presence after refining external-domain meshes.
 * @details
 * Rationale: refinement must not silently drop regions or boundaries in multi-material meshes.
 * Method: finalize external mesh, check conforming status, growth in `NE`, and nonzero counts for expected IDs.
 * If this fails: inspect `Finalize` orientation/refinement calls and any attribute mutation side effects.
 */
TEST_F(stroidTest, Finalize_ExternalDomainConformingAndRefined) {
    const auto cfg_ptr = LoadConfigFromRepo("configs/test_external_domain.toml");
    const auto& cfg = *cfg_ptr;

    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    const int initial_elements = mesh->GetNE();

    stroid::topology::Finalize(*mesh, cfg);

    EXPECT_TRUE(mesh->Conforming());
    EXPECT_GT(mesh->GetNE(), initial_elements);

    const auto volume_attr_counts = CountVolumeAttributes(*mesh);
    EXPECT_GT(volume_attr_counts.at(static_cast<int>(cfg->core_id.value())), 0);
    EXPECT_GT(volume_attr_counts.at(static_cast<int>(cfg->envelope_id.value())), 0);
    EXPECT_GT(volume_attr_counts.at(static_cast<int>(cfg->vacuum_id.value())), 0);

    const auto boundary_attr_counts = CountBoundaryAttributes(*mesh);
    EXPECT_GT(boundary_attr_counts.at(static_cast<int>(cfg->surface_bdr_id.value())), 0);
    EXPECT_GT(boundary_attr_counts.at(static_cast<int>(cfg->inf_bdr_id.value())), 0);
}

/**
 * @brief Enforces strict attribute set invariants after finalize.
 * @details
 * Rationale: presence checks alone can miss rogue IDs introduced by buggy attribute rewrites.
 * Method: assert post-finalize volume/boundary attribute keys exactly match expected sets.
 * If this fails: inspect all calls to `SetAttribute` / `SetBdrAttribute` in topology + utility code,
 * notably `src/lib/topology/topology.cpp` and `src/lib/utils/mesh_utils.cpp`.
 */
TEST_F(stroidTest, Finalize_ExternalDomainKeepsOnlyExpectedMaterialAndBoundaryIDs) {
    const auto cfg_ptr = LoadConfigFromRepo("configs/test_external_domain.toml");
    const auto& cfg = *cfg_ptr;

    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);

    stroid::topology::Finalize(*mesh, cfg);

    const auto volume_attr_counts = CountVolumeAttributes(*mesh);
    const std::set<int> expected_volume_ids = {
        static_cast<int>(cfg->core_id.value()),
        static_cast<int>(cfg->envelope_id.value()),
        static_cast<int>(cfg->vacuum_id.value())
    };
    for (const auto& [attr, count] : volume_attr_counts) {
        EXPECT_TRUE(expected_volume_ids.contains(attr));
        EXPECT_GT(count, 0);
    }
    EXPECT_EQ(volume_attr_counts.size(), expected_volume_ids.size());

    const auto boundary_attr_counts = CountBoundaryAttributes(*mesh);
    const std::set<int> expected_boundary_ids = {
        static_cast<int>(cfg->surface_bdr_id.value()),
        static_cast<int>(cfg->inf_bdr_id.value())
    };
    for (const auto& [attr, count] : boundary_attr_counts) {
        EXPECT_TRUE(expected_boundary_ids.contains(attr));
        EXPECT_GT(count, 0);
    }
    EXPECT_EQ(boundary_attr_counts.size(), expected_boundary_ids.size());
}

/**
 * @brief Verifies high-order promotion actually attaches nodal data.
 * @details
 * Rationale: projection and most quality metrics are node-based; missing nodes means pipeline misuse.
 * Method: finalize then promote, assert nodes exist and are finite.
 * If this fails: inspect `stroid::topology::PromoteToHighOrder` in
 * `src/lib/topology/curvilinear.cpp` (`H1_FECollection`, `SetNodalFESpace`).
 */
TEST_F(stroidTest, PromoteToHighOrder_SetsNodes) {
    const Config cfg;
    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    stroid::topology::Finalize(*mesh, cfg);

    stroid::topology::PromoteToHighOrder(*mesh, cfg);

    EXPECT_NE(mesh->GetNodes(), nullptr);
    EXPECT_TRUE(IsFiniteMeshNodes(*mesh));
}

/**
 * @brief Confirms projection does not produce NaN/Inf nodal coordinates.
 * @details
 * Rationale: finite nodes are the minimum numerical sanity condition for any downstream solver.
 * Method: run full pre-projection pipeline, project once, assert all node components are finite.
 * If this fails: inspect `ProjectMesh` and mapping functions in
 * `src/lib/topology/curvilinear.cpp` and `src/lib/topology/mapping.cpp`.
 */
TEST_F(stroidTest, ProjectMesh_ProducesFiniteNodes) {
    const Config cfg;
    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    stroid::topology::Finalize(*mesh, cfg);
    stroid::topology::PromoteToHighOrder(*mesh, cfg);

    stroid::topology::ProjectMesh(*mesh, cfg);

    EXPECT_TRUE(IsFiniteMeshNodes(*mesh));
}

/**
 * @brief Unit-checks equiangular mapping against closed-form tangent expression.
 * @details
 * Rationale: this catches formula or branch edits before they propagate into mesh-wide projection.
 * Method: transform a known point and compare against explicit `tan(pi/4 * ratio)` expectations.
 * If this fails: inspect `ApplyEquiangular` in `src/lib/topology/mapping.cpp`, especially dominant-axis logic.
 */
TEST_F(stroidTest, ApplyEquiangular_BasicTransform) {
    mfem::Vector pos(3);
    pos(0) = 1.0;
    pos(1) = 0.5;
    pos(2) = -0.25;

    stroid::topology::ApplyEquiangular(pos);

    const double expected_y = 1.0 * std::tan(kPi / 4.0 * (0.5 / 1.0));
    const double expected_z = 1.0 * std::tan(kPi / 4.0 * (-0.25 / 1.0));
    EXPECT_NEAR(pos(1), expected_y, 1e-12);
    EXPECT_NEAR(pos(2), expected_z, 1e-12);
}

/**
 * @brief Verifies spheroidal flattening scales only the z component as configured.
 * @details
 * Rationale: flattening is intentionally simple and should remain easy to reason about.
 * Method: load flattening config, apply transform to z-axis point, check exact expected z.
 * If this fails: inspect `ApplySpheroidal` in `src/lib/topology/mapping.cpp` and fixture values in
 * `configs/test_flattening.toml`.
 */
TEST_F(stroidTest, ApplySpheroidal_FlattensZ) {
    const auto cfg_ptr = LoadConfigFromRepo("configs/test_flattening.toml");
    const auto& cfg = *cfg_ptr;

    mfem::Vector pos(3);
    pos(0) = 0.0;
    pos(1) = 0.0;
    pos(2) = 10.0;

    stroid::topology::ApplySpheroidal(pos, cfg);

    EXPECT_NEAR(pos(2), 8.0, 1e-12);
}

/**
 * @brief Ensures axis-aligned points inside the core stay unchanged in this mapping regime.
 * @details
 * Rationale: axis points are symmetry anchors; any drift here is usually a serious branch regression.
 * Method: map `(1,0,0)` with default config and assert identity.
 * If this fails: inspect `TransformPoint` core-zone blend path and instability guard in
 * `src/lib/topology/mapping.cpp`.
 */
TEST_F(stroidTest, TransformPoint_AxisInsideCore_NoChange) {
    const Config cfg;

    mfem::Vector pos(3);
    pos(0) = 1.0;
    pos(1) = 0.0;
    pos(2) = 0.0;

    stroid::topology::TransformPoint(pos, cfg, 0);

    EXPECT_NEAR(pos(0), 1.0, 1e-12);
    EXPECT_NEAR(pos(1), 0.0, 1e-12);
    EXPECT_NEAR(pos(2), 0.0, 1e-12);
}

/**
 * @brief Ensures axis-aligned envelope points remain unchanged under default spherical setup.
 * @details
 * Rationale: this verifies envelope branch consistency and avoids silent radial drift.
 * Method: map `(3,0,0)` and assert identity under default parameters.
 * If this fails: inspect `TransformPoint` envelope branch and normalized-direction reconstruction logic.
 */
TEST_F(stroidTest, TransformPoint_AxisEnvelope_NoChange) {
    const Config cfg;

    mfem::Vector pos(3);
    pos(0) = 3.0;
    pos(1) = 0.0;
    pos(2) = 0.0;

    stroid::topology::TransformPoint(pos, cfg, 0);

    EXPECT_NEAR(pos(0), 3.0, 1e-12);
    EXPECT_NEAR(pos(1), 0.0, 1e-12);
    EXPECT_NEAR(pos(2), 0.0, 1e-12);
}

/**
 * @brief Checks continuity near `r_core` and `r_star` transition surfaces.
 * @details
 * Rationale: discontinuities at interfaces destabilize high-order interpolation and integration.
 * Method: evaluate points at `r*(1±eps)` and assert mapped separation stays small in L2 norm.
 * If this fails: inspect transition formulas in `TransformPoint` (core blend, envelope mapping)
 * and any recent edits to `core_steepness` handling.
 */
TEST_F(stroidTest, TransformPoint_IsContinuousAcrossCoreAndStarInterfaces) {
    const Config cfg;
    constexpr double eps = 1e-6;

    mfem::Vector dir(3);
    dir(0) = 1.0;
    dir(1) = 0.6;
    dir(2) = -0.4;

    mfem::Vector near_core_left = dir;
    near_core_left *= cfg->r_core.value() * (1.0 - eps);
    mfem::Vector near_core_right = dir;
    near_core_right *= cfg->r_core.value() * (1.0 + eps);

    const mfem::Vector core_left_mapped = TransformCopy(near_core_left, cfg);
    const mfem::Vector core_right_mapped = TransformCopy(near_core_right, cfg);

    mfem::Vector diff = core_left_mapped;
    diff -= core_right_mapped;
    EXPECT_LT(diff.Norml2(), 1e-3);

    mfem::Vector near_star_left = dir;
    near_star_left *= cfg->r_star.value() * (1.0 - eps);
    mfem::Vector near_star_right = dir;
    near_star_right *= cfg->r_star.value() * (1.0 + eps);

    const mfem::Vector star_left_mapped = TransformCopy(near_star_left, cfg);
    const mfem::Vector star_right_mapped = TransformCopy(near_star_right, cfg);

    diff = star_left_mapped;
    diff -= star_right_mapped;
    EXPECT_LT(diff.Norml2(), 1e-3);
}

/**
 * @brief Verifies expected sign and axis-permutation symmetry in the spherical case.
 * @details
 * Rationale: symmetry violations usually indicate branch asymmetry bugs in mapping logic.
 * Method: compare mapped values for original, sign-flipped, and axis-swapped points.
 * If this fails: inspect dominant-axis branching in `ApplyEquiangular` and normalization flow in
 * `TransformPoint`.
 */
TEST_F(stroidTest, TransformPoint_RespectsSignAndAxisPermutationSymmetryWithoutFlattening) {
    const Config cfg;

    mfem::Vector p(3);
    p(0) = 4.0;
    p(1) = 2.0;
    p(2) = 1.0;

    mfem::Vector p_neg = p;
    p_neg *= -1.0;

    mfem::Vector p_swapped(3);
    p_swapped(0) = p(1);
    p_swapped(1) = p(0);
    p_swapped(2) = p(2);

    const mfem::Vector mapped = TransformCopy(p, cfg);
    const mfem::Vector mapped_neg = TransformCopy(p_neg, cfg);
    const mfem::Vector mapped_swapped = TransformCopy(p_swapped, cfg);

    EXPECT_NEAR(mapped_neg(0), -mapped(0), 1e-12);
    EXPECT_NEAR(mapped_neg(1), -mapped(1), 1e-12);
    EXPECT_NEAR(mapped_neg(2), -mapped(2), 1e-12);

    EXPECT_NEAR(mapped_swapped(0), mapped(1), 1e-12);
    EXPECT_NEAR(mapped_swapped(1), mapped(0), 1e-12);
    EXPECT_NEAR(mapped_swapped(2), mapped(2), 1e-12);
}

/**
 * @brief Smoke-tests mesh serialization to MFEM format.
 * @details
 * Rationale: I/O regressions are easy to miss during geometry-focused development.
 * Method: write a finalized mesh to temp storage and assert file exists and is non-empty.
 * If this fails: inspect `stroid::IO::SaveMesh` in `src/lib/IO/mesh.cpp` and local filesystem perms.
 */
TEST_F(stroidTest, SaveMesh_WritesFile) {
    const Config cfg;
    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    stroid::topology::Finalize(*mesh, cfg);

    const std::filesystem::path tmp_dir = std::filesystem::temp_directory_path();
    const std::filesystem::path mesh_path = tmp_dir / "stroid_test_mesh.mesh";

    stroid::IO::SaveMesh(*mesh, mesh_path.string());

    ASSERT_TRUE(std::filesystem::exists(mesh_path));
    EXPECT_GT(std::filesystem::file_size(mesh_path), 0u);

    std::error_code ec;
    std::filesystem::remove(mesh_path, ec);
}

/**
 * @brief End-to-end baseline pipeline test for the default domain.
 * @details
 * Rationale: validates the canonical operation order used by both library examples and CLI.
 * Method: execute full pipeline and assert non-empty, nodal, finite output mesh.
 * If this fails: check call-order assumptions and recent edits in
 * `src/lib/topology/topology.cpp` / `src/lib/topology/curvilinear.cpp`.
 */
TEST_F(stroidTest, EndToEnd_BuildFinalizePromoteProject) {
    const Config cfg;
    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    stroid::topology::Finalize(*mesh, cfg);
    stroid::topology::PromoteToHighOrder(*mesh, cfg);
    stroid::topology::ProjectMesh(*mesh, cfg);

    EXPECT_GT(mesh->GetNE(), 0);
    EXPECT_NE(mesh->GetNodes(), nullptr);
    EXPECT_TRUE(IsFiniteMeshNodes(*mesh));
}

/**
 * @brief End-to-end pipeline test for external-domain meshes.
 * @details
 * Rationale: confirms the same pipeline remains valid when vacuum blocks are included.
 * Method: run full external-domain pipeline and assert conforming + finite nodal output.
 * If this fails: inspect external-domain topology assembly and projection loops over mixed attributes.
 */
TEST_F(stroidTest, EndToEnd_ExternalDomainBuildFinalizePromoteProject) {
    const auto cfg_ptr= LoadConfigFromRepo("configs/test_external_domain.toml");
    const auto& cfg = *cfg_ptr;

    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    stroid::topology::Finalize(*mesh, cfg);
    stroid::topology::PromoteToHighOrder(*mesh, cfg);
    stroid::topology::ProjectMesh(*mesh, cfg);

    EXPECT_TRUE(mesh->Conforming());
    EXPECT_GT(mesh->GetNE(), 0);
    EXPECT_NE(mesh->GetNodes(), nullptr);
    EXPECT_TRUE(IsFiniteMeshNodes(*mesh));
}

/**
 * @brief Verifies stellar volume invariance with respect to adding a vacuum shell.
 * @details
 * Rationale: the vacuum region should extend domain extent, not alter stellar mass volume.
 * Method: integrate core+envelope volume in both configs and compare relative difference.
 * If this fails: inspect attribute-filtered volume helpers in this file and mapping/topology changes
 * that may leak starside nodes into vacuum geometry.
 */
TEST_F(stroidTest, Volume_StellarDomainMatchesWithAndWithoutExternalDomain) {
    const auto no_external_cfg_ptr = LoadConfigFromRepo("configs/test_volume_no_external.toml");
    const auto with_external_cfg_ptr = LoadConfigFromRepo("configs/test_volume_with_external.toml");

    const auto& no_external_cfg = *no_external_cfg_ptr;
    const auto& with_external_cfg = *with_external_cfg_ptr;


    const std::unique_ptr<mfem::Mesh> no_external_mesh = stroid::topology::BuildSkeleton(no_external_cfg);
    stroid::topology::Finalize(*no_external_mesh, no_external_cfg);
    stroid::topology::PromoteToHighOrder(*no_external_mesh, no_external_cfg);
    stroid::topology::ProjectMesh(*no_external_mesh, no_external_cfg);

    const std::unique_ptr<mfem::Mesh> with_external_mesh = stroid::topology::BuildSkeleton(with_external_cfg);
    stroid::topology::Finalize(*with_external_mesh, with_external_cfg);
    stroid::topology::PromoteToHighOrder(*with_external_mesh, with_external_cfg);
    stroid::topology::ProjectMesh(*with_external_mesh, with_external_cfg);

    const std::set<int> stellar_attrs_no_external = {
        static_cast<int>(no_external_cfg->core_id.value()),
        static_cast<int>(no_external_cfg->envelope_id.value())
    };
    const std::set<int> stellar_attrs_with_external = {
        static_cast<int>(with_external_cfg->core_id.value()),
        static_cast<int>(with_external_cfg->envelope_id.value())
    };

    const double stellar_volume_no_external = ComputeMeshVolumeForAttributes(*no_external_mesh, stellar_attrs_no_external);
    const double stellar_volume_with_external = ComputeMeshVolumeForAttributes(*with_external_mesh, stellar_attrs_with_external);

    const double rel_diff = std::abs(stellar_volume_with_external - stellar_volume_no_external) /
                            std::max(stellar_volume_with_external, stellar_volume_no_external);
    EXPECT_LT(rel_diff, 5e-3);
}

/**
 * @brief Confirms total volume decomposition into stellar + vacuum components.
 * @details
 * Rationale: this explicitly checks that vacuum exclusion logic is doing real work, not a no-op.
 * Method: on external mesh, compute total, stellar-only, and vacuum-only volumes and enforce
 * additive consistency.
 * If this fails: inspect `ComputeMeshVolume*` helpers and region attribute IDs in config fixtures.
 */
TEST_F(stroidTest, Volume_ExternalMeshExcludesVacuumWhenRequested) {
    const auto cfg_ptr = LoadConfigFromRepo("configs/test_volume_with_external.toml");
    const auto& cfg = *cfg_ptr;


    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    stroid::topology::Finalize(*mesh, cfg);
    stroid::topology::PromoteToHighOrder(*mesh, cfg);
    stroid::topology::ProjectMesh(*mesh, cfg);

    const std::set<int> stellar_attrs = {
        static_cast<int>(cfg->core_id.value()),
        static_cast<int>(cfg->envelope_id.value())
    };
    const std::set<int> vacuum_attr = {static_cast<int>(cfg->vacuum_id.value())};

    const double total_volume = ComputeMeshVolume(*mesh);
    const double stellar_volume = ComputeMeshVolumeForAttributes(*mesh, stellar_attrs);
    const double vacuum_volume = ComputeMeshVolumeForAttributes(*mesh, vacuum_attr);

    EXPECT_GT(vacuum_volume, 0.0);
    EXPECT_GT(total_volume, stellar_volume);
    EXPECT_NEAR(total_volume, stellar_volume + vacuum_volume, total_volume * 1e-9 + 1e-12);
}

/**
 * @brief Compares direct Jacobian-based stellar volume to analytic sphere volume.
 * @details
 * Rationale: anchors numerical integration against a closed-form reference in the spherical limit.
 * Method: integrate core+envelope using element transformations, compare to `4/3*pi*r_star^3`.
 * If this fails: inspect mapping spherical path (`flattening=0`) and quadrature order in
 * `IntegrateElementVolume`.
 */
TEST_F(stroidTest, Volume_SphericalStellarDomainMatchesAnalyticSphere) {
    const auto cfg_ptr = LoadConfigFromRepo("configs/test_volume_spherical_no_external.toml");
    const auto& cfg = *cfg_ptr;

    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    stroid::topology::Finalize(*mesh, cfg);
    stroid::topology::PromoteToHighOrder(*mesh, cfg);
    stroid::topology::ProjectMesh(*mesh, cfg);

    const std::set<int> stellar_attrs = {
        static_cast<int>(cfg->core_id.value()),
        static_cast<int>(cfg->envelope_id.value())
    };

    const double measured_volume = ComputeMeshVolumeForAttributes(*mesh, stellar_attrs);
    const double analytic_volume = 4.0 / 3.0 * kPi * std::pow(cfg->r_star.value(), 3.0);
    const double rel_err = std::abs(measured_volume - analytic_volume) / analytic_volume;

    EXPECT_LT(rel_err, 1e-2);
}

/**
 * @brief Repeats spherical analytic-volume check via MFEM `DomainLFIntegrator`.
 * @details
 * Rationale: independent integration machinery lowers the risk of helper-specific false confidence.
 * Method: build attribute-weighted `PWConstCoefficient` (core+envelope=1, vacuum=0), assemble
 * domain linear form, and compare against analytic sphere volume.
 * If this fails: inspect coefficient indexing (attr-1 convention), `ComputeStellarVolumeWithDomainLFIntegrator`,
 * and MFEM assembly setup in this test file.
 */
TEST_F(stroidTest, Volume_SphericalStellarDomainDomainLFIntegratorMatchesAnalyticSphere) {
    const auto cfg_ptr = LoadConfigFromRepo("configs/test_volume_spherical_with_external.toml");
    const auto& cfg = *cfg_ptr;

    std::unique_ptr<mfem::Mesh> mesh = BuildProjectedMesh(cfg);
    const double measured_volume = ComputeStellarVolumeWithDomainLFIntegrator(*mesh, cfg);
    const double analytic_volume = 4.0 / 3.0 * kPi * std::pow(cfg->r_star.value(), 3.0);
    const double rel_err = std::abs(measured_volume - analytic_volume) / analytic_volume;

    EXPECT_LT(rel_err, 1e-2);
}

/**
 * @brief Enforces baseline conditioning bounds for the default projected mesh.
 * @details
 * Rationale: this guards against silent degradation in element quality that may still pass finiteness checks.
 * Method: sample Jacobian stats over quadrature points and assert positivity + distortion/stretch bounds.
 * If this fails: inspect mapping formulas in `src/lib/topology/mapping.cpp` and any changes to
 * refinement/order config used by `configs/test_volume_spherical_no_external.toml`.
 */
TEST_F(stroidTest, Conditioning_DefaultMeshHasPositiveJacobiansAndReasonableShape) {
    const auto cfg_ptr = LoadConfigFromRepo("configs/test_volume_spherical_no_external.toml");
    const auto& cfg = *cfg_ptr;

    const std::unique_ptr<mfem::Mesh> mesh = BuildProjectedMesh(cfg);

    const ConditioningStats stats = CollectConditioningStats(*mesh, {});

    ASSERT_GT(stats.samples, 0);
    EXPECT_GT(stats.min_det, 1e-10);
    EXPECT_LT(stats.max_det / stats.min_det, 1e6);
    EXPECT_GT(stats.min_scaled_jac, 1e-3);
    EXPECT_LT(stats.max_stretch_ratio, 50.0);
    EXPECT_LT(stats.max_edge_ratio, 50.0);
}

/**
 * @brief Applies Jacobian conditioning checks independently to core, envelope, and vacuum regions.
 * @details
 * Rationale: global stats can hide localized failures; region-level checks make regressions diagnosable.
 * Method: collect conditioning statistics per attribute and enforce positive Jacobians + scaled-Jacobian floors.
 * If this fails: inspect region-specific mapping behavior in `TransformPoint` and verify attribute
 * assignment in `BuildSkeleton`.
 */
TEST_F(stroidTest, Conditioning_ExternalMeshPerRegionHasPositiveJacobians) {
    const auto cfg_ptr = LoadConfigFromRepo("configs/test_volume_spherical_with_external.toml");
    const auto& cfg = *cfg_ptr;

    const std::unique_ptr<mfem::Mesh> mesh = BuildProjectedMesh(cfg);

    const ConditioningStats core_stats = CollectConditioningStats(*mesh, {static_cast<int>(cfg->core_id.value())});
    const ConditioningStats envelope_stats = CollectConditioningStats(*mesh, {static_cast<int>(cfg->envelope_id.value())});
    const ConditioningStats vacuum_stats = CollectConditioningStats(*mesh, {static_cast<int>(cfg->vacuum_id.value())});

    ASSERT_GT(core_stats.samples, 0);
    ASSERT_GT(envelope_stats.samples, 0);
    ASSERT_GT(vacuum_stats.samples, 0);

    EXPECT_GT(core_stats.min_det, 1e-10);
    EXPECT_GT(envelope_stats.min_det, 1e-10);
    EXPECT_GT(vacuum_stats.min_det, 1e-10);

    EXPECT_GT(core_stats.min_scaled_jac, 1e-3);
    EXPECT_GT(envelope_stats.min_scaled_jac, 2e-2);
    EXPECT_GT(vacuum_stats.min_scaled_jac, 1e-3);
}

/**
 * @brief Validates orientation quality via flipped-element and flipped-boundary markers.
 * @details
 * Rationale: this is a direct orientation sanity check using project utilities already used for debugging.
 * Method: run `MarkFlippedElements`/`MarkFlippedBoundaryElements` and assert sentinel attrs are absent.
 * If this fails: inspect Jacobian sign behavior and boundary normal orientation code in
 * `src/lib/utils/mesh_utils.cpp`, then trace upstream mapping changes.
 */
TEST_F(stroidTest, Conditioning_DefaultMeshHasNoFlippedElementsOrBoundaryFaces) {
    const auto cfg_ptr = LoadConfigFromRepo("configs/test_volume_spherical_no_external.toml");
    const auto& cfg = *cfg_ptr;

    std::unique_ptr<mfem::Mesh> mesh = BuildProjectedMesh(cfg);

    stroid::utils::MarkFlippedElements(*mesh);
    stroid::utils::MarkFlippedBoundaryElements(*mesh);

    const auto volume_attr_counts = CountVolumeAttributes(*mesh);
    const auto boundary_attr_counts = CountBoundaryAttributes(*mesh);

    EXPECT_FALSE(volume_attr_counts.contains(999));
    EXPECT_FALSE(boundary_attr_counts.contains(500));
}

TEST_F(stroidTest, PolynomainalProjection) {
    const auto cfg_ptr = LoadConfigFromRepo("configs/test_polynomial_projection.toml");
    const auto& cfg = *cfg_ptr;

    std::unique_ptr<mfem::Mesh> mesh = BuildProjectedMesh(cfg);

    const int geom_order = mesh->GetNodes()->FESpace()->GetMaxElementOrder();
    const int space_dim = mesh->Dimension();

    mfem::H1_FECollection fec(geom_order, space_dim);
    mfem::FiniteElementSpace fes(mesh.get(), &fec);

    auto ProjectedFunction = [](const mfem::Vector& x) {
        const double r = x.Norml2();
        return 1 + 7 * r * r - 2 * r;
    };

    mfem::GridFunction projected_u(&fes);
    mfem::FunctionCoefficient u_coeff(ProjectedFunction);
    projected_u.ProjectCoefficient(u_coeff);

    mfem::Vector x(space_dim);
    x = 0.0;
    for (double t = 0; t <= 1; t+= 0.01) {
        x(0) = t;

        double analytic_val = ProjectedFunction(x);
        double projected_val = EvalGridFunctionAtPoint(*mesh, x, projected_u).value_or(std::numeric_limits<double>::quiet_NaN());

        double rel_err = std::abs(projected_val - analytic_val) / analytic_val;
        EXPECT_LT(rel_err, 1e-12);
    }
}

TEST_F(stroidTest, TranscendtalProjection) {
    const auto cfg_ptr = LoadConfigFromRepo("configs/test_polynomial_projection.toml");
    const auto& cfg = *cfg_ptr;

    std::unique_ptr<mfem::Mesh> mesh = BuildProjectedMesh(cfg);

    const int geom_order = mesh->GetNodes()->FESpace()->GetMaxElementOrder();
    const int space_dim = mesh->Dimension();

    mfem::H1_FECollection fec(geom_order, space_dim);
    mfem::FiniteElementSpace fes(mesh.get(), &fec);


    auto ProjectedFunction = [](const mfem::Vector& x) {
        const double r = x.Norml2();
        if (r <= 1e-8) return 1.0;
        return std::sin(r)/r;
    };

    auto expansion = [](const double t, const int order) {
        double val = 0.0;
        for (int k = 0; k < order; ++k) {
            const double sign = (k % 2 == 0) ? 1.0 : -1.0;
            const double term = sign * std::pow(t, 2 * k) / std::tgamma(2 * k + 2);
            val += term;
        }
        return val;
    };

    auto expansion_err = [geom_order, &expansion](const double r) {
        const double expansion_val = expansion(r, geom_order);

        const double analytic_val = std::sin(r)/r;
        return std::abs((expansion_val - analytic_val))/std::abs(analytic_val);
    };

    double max_estimated_truncation_error = 0.0;
    for (double t = 0; t < 1; t+= 0.01) {
        double trunc_err = expansion_err(t);
        max_estimated_truncation_error = std::max(max_estimated_truncation_error, trunc_err);
    }

    mfem::GridFunction projected_u(&fes);
    mfem::FunctionCoefficient u_coeff(ProjectedFunction);
    projected_u.ProjectCoefficient(u_coeff);

    mfem::Vector x(space_dim);
    x = 0.0;
    for (double t = 0; t <= 1; t+= 0.01) {
        x(0) = t;

        double analytic_val = ProjectedFunction(x);
        double projected_val = EvalGridFunctionAtPoint(*mesh, x, projected_u).value_or(std::numeric_limits<double>::quiet_NaN());

        double rel_err = std::abs(projected_val - analytic_val) / analytic_val;
        EXPECT_LT(rel_err, 10*max_estimated_truncation_error);
    }
}

TEST_F(stroidTest, Refinement_UniformRefinementProducesExpectedElementCounts) {
    const auto cfg_ptr = LoadConfigFromRepo("configs/test_volume_spherical_no_external.toml");
    const auto& cfg = *cfg_ptr;

    stroid::StroidMesh mesh;
    EXPECT_NO_THROW(mesh = stroid::GenerateMesh(cfg));
    size_t init_elements = mesh.mesh->GetNE();

    stroid::refinement::UniformRefinement(mesh, 1);
    EXPECT_EQ(mesh.mesh->GetNE(), init_elements * 8);
}

TEST_F(stroidTest, ExteriorCoordinate_HasValidRangeAndExactBoundaryTraces) {
    const auto cfg_ptr = LoadConfigFromRepo("configs/test_volume_with_external.toml");
    const auto& cfg = *cfg_ptr;

    stroid::StroidMesh mesh;
    ASSERT_NO_THROW(mesh = stroid::GenerateMesh(cfg));

    ExpectExteriorCoordinateRange(mesh);
    ExpectExteriorCoordinateBoundaryTraces(mesh);
}

TEST_F(stroidTest, ExteriorCoordinate_IsRebuiltAfterUniformRefinement) {
    const auto cfg_ptr = LoadConfigFromRepo("configs/test_volume_with_external.toml");
    const auto& cfg = *cfg_ptr;

    stroid::StroidMesh mesh;
    ASSERT_NO_THROW(mesh = stroid::GenerateMesh(cfg));
    ASSERT_NE(mesh.exterior_coordinate, nullptr);

    const int initial_elements = mesh.mesh->GetNE();
    const int initial_coordinate_dofs = mesh.exterior_coordinate->space->GetNDofs();

    ASSERT_NO_THROW(stroid::refinement::UniformRefinement(mesh, 1));
    ASSERT_NE(mesh.exterior_coordinate, nullptr);
    EXPECT_EQ(mesh.mesh->GetNE(), initial_elements * 8);
    EXPECT_GT(mesh.exterior_coordinate->space->GetNDofs(), initial_coordinate_dofs);

    ExpectExteriorCoordinateRange(mesh);
    ExpectExteriorCoordinateBoundaryTraces(mesh);
}

TEST_F(stroidTest, ExteriorCoordinate_SurvivesSaveAndLoad) {
    const auto cfg_ptr = LoadConfigFromRepo("configs/test_volume_with_external.toml");
    const auto& cfg = *cfg_ptr;

    stroid::StroidMesh original;
    ASSERT_NO_THROW(original = stroid::GenerateMesh(cfg));
    ASSERT_NE(original.exterior_coordinate, nullptr);

    const std::filesystem::path path = std::filesystem::temp_directory_path() / "stroid_exterior_coordinate_round_trip.smesh";
    ASSERT_NO_THROW(stroid::IO::SaveStroidMesh(original, path.string(), "Exterior-coordinate round-trip test"));

    auto loaded_result = stroid::IO::LoadStroidMesh(path.string());
    if (!loaded_result.has_value()) FAIL() << loaded_result.error();
    stroid::StroidMesh loaded = std::move(*loaded_result);

    ASSERT_NE(loaded.exterior_coordinate, nullptr);
    ASSERT_EQ(loaded.exterior_coordinate->space->GetNDofs(), original.exterior_coordinate->space->GetNDofs());
    ASSERT_EQ(loaded.exterior_coordinate->values->Size(), original.exterior_coordinate->values->Size());

    for (int dof = 0; dof < original.exterior_coordinate->values->Size(); ++dof) {
        EXPECT_DOUBLE_EQ((*loaded.exterior_coordinate->values)(dof), (*original.exterior_coordinate->values)(dof));
    }

    ExpectExteriorCoordinateRange(loaded);
    ExpectExteriorCoordinateBoundaryTraces(loaded);

    std::error_code error;
    std::filesystem::remove(path, error);
    EXPECT_FALSE(error);
}

TEST_F(stroidTest, ExteriorCoordinate_IsAbsentWithoutExternalDomainAcrossSaveAndLoad) {
    const auto cfg_ptr = LoadConfigFromRepo("configs/test_volume_spherical_no_external.toml");
    const auto& cfg = *cfg_ptr;

    stroid::StroidMesh original;
    ASSERT_NO_THROW(original = stroid::GenerateMesh(cfg));
    EXPECT_EQ(original.exterior_coordinate, nullptr);

    const std::filesystem::path path = std::filesystem::temp_directory_path() / "stroid_no_exterior_coordinate_round_trip.smesh";
    ASSERT_NO_THROW(stroid::IO::SaveStroidMesh(original, path.string(), "No-exterior-coordinate round-trip test"));

    auto loaded_result = stroid::IO::LoadStroidMesh(path.string());
    if (!loaded_result.has_value()) FAIL() << loaded_result.error();
    EXPECT_EQ(loaded_result->exterior_coordinate, nullptr);

    std::error_code error;
    std::filesystem::remove(path, error);
    EXPECT_FALSE(error);
}

TEST_F(stroidTest, ExteriorCoordinate_IsReconstructedWhenLoadingLegacyFiles) {
    const auto cfg_ptr = LoadConfigFromRepo("configs/test_volume_with_external.toml");
    const auto& cfg = *cfg_ptr;

    stroid::StroidMesh original;
    ASSERT_NO_THROW(original = stroid::GenerateMesh(cfg));
    ASSERT_NE(original.exterior_coordinate, nullptr);

    const std::filesystem::path path = std::filesystem::temp_directory_path() / "stroid_legacy_exterior_coordinate.smesh";
    ASSERT_NO_THROW(stroid::IO::SaveStroidMesh(original, path.string(), "Legacy exterior-coordinate reconstruction test"));

    std::ifstream input(path);
    ASSERT_TRUE(input.is_open());
    std::string contents((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());

    constexpr std::string_view begin_marker = "BEGIN BLOCK EXTERIOR_COORDINATE";
    constexpr std::string_view end_marker = "END BLOCK EXTERIOR_COORDINATE";
    const size_t begin = contents.find(begin_marker);
    const size_t end_begin = contents.find(end_marker);
    ASSERT_NE(begin, std::string::npos);
    ASSERT_NE(end_begin, std::string::npos);

    size_t end = end_begin + end_marker.size();
    if (end < contents.size() && contents[end] == '\n') ++end;
    contents.erase(begin, end - begin);

    std::istringstream legacy_stream(contents);
    auto loaded_result = stroid::IO::ParseStroidMesh(legacy_stream);
    if (!loaded_result.has_value()) FAIL() << loaded_result.error();

    stroid::StroidMesh loaded = std::move(*loaded_result);
    ASSERT_NE(loaded.exterior_coordinate, nullptr);
    ASSERT_EQ(loaded.exterior_coordinate->values->Size(), original.exterior_coordinate->values->Size());

    for (int dof = 0; dof < original.exterior_coordinate->values->Size(); ++dof) {
        EXPECT_NEAR((*loaded.exterior_coordinate->values)(dof), (*original.exterior_coordinate->values)(dof), 1.0e-12);
    }

    ExpectExteriorCoordinateRange(loaded);
    ExpectExteriorCoordinateBoundaryTraces(loaded);

    std::error_code error;
    std::filesystem::remove(path, error);
    EXPECT_FALSE(error);
}

TEST_F(stroidTest, Stats_ComputeStats) {
    const auto cfg_ptr = LoadConfigFromRepo("configs/test_volume_with_external.toml");
    const auto& cfg = *cfg_ptr;

    stroid::StroidMesh mesh;
    EXPECT_NO_THROW(mesh = stroid::GenerateMesh(cfg));

    stroid::stats::MeshStats stats = stroid::stats::ComputeMeshStats(mesh);
    std::println("{}", stats);

}

namespace {

std::unique_ptr<Config> MultiBlockConfiguration(int order, int refinement, bool external, double flattening = 0.0) {
    auto cfg = std::make_unique<Config>();
    cfg->mutate([&](stroid::config::MeshConfig& value) {
        value.core_mapping = "multi_block";
        value.order = order;
        value.refinement_levels = refinement;
        value.include_external_domain = external;
        value.flattening = flattening;
        value.optimization_methods = stroid::config::OptimizationMethods{false, false};
    });
    return cfg;
}

// Unlike CollectConditioningStats, this uses the signed determinant, actual
// singular values, and a closed sample grid including vertices/edges/faces.
// Column-length ratios and open quadrature points miss the old core-corner defect.
void ExpectClosedGridCoreConditioning(mfem::Mesh& mesh, int coreAttribute, double maximumCondition = 10.0) {
    int coreElements = 0;
    double largestCondition = 0.0;
    double smallestDeterminant = std::numeric_limits<double>::infinity();
    for (int element = 0; element < mesh.GetNE(); ++element) {
        if (mesh.GetAttribute(element) != coreAttribute) continue;
        ++coreElements;
        auto* transformation = mesh.GetElementTransformation(element);
        ASSERT_EQ(transformation->GetGeometryType(), mfem::Geometry::CUBE);
        for (double x : {0.0, 0.01, 0.5, 0.99, 1.0}) {
            for (double y : {0.0, 0.01, 0.5, 0.99, 1.0}) {
                for (double z : {0.0, 0.01, 0.5, 0.99, 1.0}) {
                    mfem::IntegrationPoint point;
                    point.Set3(x, y, z);
                    transformation->SetIntPoint(&point);
                    const auto& jacobian = transformation->Jacobian();
                    const double determinant = jacobian.Det();
                    const double minimumSingular = jacobian.CalcSingularvalue(2);
                    const double maximumSingular = jacobian.CalcSingularvalue(0);
                    ASSERT_TRUE(std::isfinite(determinant));
                    ASSERT_GT(determinant, 0.0) << "element=" << element << " point=" << x << ',' << y << ',' << z;
                    ASSERT_TRUE(std::isfinite(minimumSingular));
                    ASSERT_GT(minimumSingular, 0.0) << "element=" << element;
                    const double condition = maximumSingular / minimumSingular;
                    ASSERT_TRUE(std::isfinite(condition));
                    ASSERT_LT(condition, maximumCondition)
                        << "element=" << element << " point=" << x << ',' << y << ',' << z;
                    smallestDeterminant = std::min(smallestDeterminant, determinant);
                    largestCondition = std::max(largestCondition, condition);
                }
            }
        }
    }
    EXPECT_GT(coreElements, 0);
    EXPECT_GT(smallestDeterminant, 0.0);
    EXPECT_LT(largestCondition, maximumCondition);
}

void ExpectCoreFaceContinuity(mfem::Mesh& mesh, int coreAttribute) {
    int faces = 0;
    mfem::Vector left(3), right(3);
    for (int face = 0; face < mesh.GetNumFaces(); ++face) {
        auto* transformation = mesh.GetFaceElementTransformations(face);
        if (transformation == nullptr || transformation->Elem1 == nullptr || transformation->Elem2 == nullptr) continue;
        if (transformation->Elem1->Attribute != coreAttribute && transformation->Elem2->Attribute != coreAttribute) continue;
        ++faces;
        for (double x : {0.0, 0.25, 0.5, 0.75, 1.0}) {
            for (double y : {0.0, 0.25, 0.5, 0.75, 1.0}) {
                mfem::IntegrationPoint point;
                point.Set2(x, y);
                transformation->SetAllIntPoints(&point);
                transformation->Elem1->Transform(transformation->Elem1->GetIntPoint(), left);
                transformation->Elem2->Transform(transformation->Elem2->GetIntPoint(), right);
                left -= right;
                EXPECT_LT(left.Norml2(), 2.0e-12) << "face=" << face;
            }
        }
    }
    EXPECT_GT(faces, 0);
}

} // namespace

TEST_F(stroidTest, MultiBlockCore_DefaultTopologyCountsAndAttributes) {
    EXPECT_EQ(stroid::config::MeshConfig{}.core_mapping.value(), "multi_block");
    for (const bool external : {false, true}) {
        SCOPED_TRACE(external);
        auto cfg = MultiBlockConfiguration(2, 0, external);
        auto mesh = stroid::topology::BuildSkeleton(*cfg);
        ASSERT_NE(mesh, nullptr);
        EXPECT_EQ(mesh->GetNV(), external ? 32 : 24);
        EXPECT_EQ(mesh->GetNE(), external ? 19 : 13);
        EXPECT_EQ(mesh->GetNBE(), external ? 12 : 6);
        const auto volumes = CountVolumeAttributes(*mesh);
        EXPECT_EQ(volumes.at(1), 7);
        EXPECT_EQ(volumes.at(2), 6);
        EXPECT_EQ(volumes.contains(3), external);
        if (external) EXPECT_EQ(volumes.at(3), 6);
        const auto boundaries = CountBoundaryAttributes(*mesh);
        EXPECT_EQ(boundaries.at(1), 6);
        EXPECT_EQ(boundaries.contains(2), external);
        if (external) EXPECT_EQ(boundaries.at(2), 6);

        cfg->mutate([](stroid::config::MeshConfig& value) { value.core_mapping = "spherified"; });
        auto legacy = stroid::topology::BuildSkeleton(*cfg);
        EXPECT_EQ(legacy->GetNE(), external ? 13 : 7);
        EXPECT_EQ(CountVolumeAttributes(*legacy).at(1), 1);
    }
}

TEST_F(stroidTest, MultiBlockCore_RejectsUnknownMappingAndInvalidGeometryConfiguration) {
    auto cfg = MultiBlockConfiguration(2, 0, true);
    cfg->mutate([](stroid::config::MeshConfig& value) { value.core_mapping = "not_a_core_mapping"; });
    EXPECT_THROW(stroid::topology::BuildSkeleton(*cfg), std::invalid_argument);
    cfg = MultiBlockConfiguration(2, 0, true);
    cfg->mutate([](stroid::config::MeshConfig& value) { value.r_core = value.r_star; });
    EXPECT_THROW(stroid::topology::BuildSkeleton(*cfg), std::invalid_argument);
    cfg = MultiBlockConfiguration(2, 0, true);
    cfg->mutate([](stroid::config::MeshConfig& value) { value.r_infinity = value.r_star; });
    EXPECT_THROW(stroid::topology::BuildSkeleton(*cfg), std::invalid_argument);
    cfg = MultiBlockConfiguration(2, 0, true);
    cfg->mutate([](stroid::config::MeshConfig& value) { value.flattening = 1.0; });
    EXPECT_THROW(stroid::topology::BuildSkeleton(*cfg), std::invalid_argument);
}

TEST_F(stroidTest, MultiBlockCore_MapHasAffineInnerCubeAndContinuousSphericalInterface) {
    auto cfg = MultiBlockConfiguration(4, 0, true);
    const double radius = (*cfg)->r_core.value();
    for (int axis = 0; axis < 3; ++axis) {
        for (double sign : {-1.0, 1.0}) {
            for (double a : {-1.0, -0.4, 0.0, 0.6, 1.0}) {
                for (double b : {-1.0, -0.3, 0.0, 0.7, 1.0}) {
                    mfem::Vector direction(3);
                    direction(axis) = sign;
                    direction((axis + 1) % 3) = a;
                    direction((axis + 2) % 3) = b;
                    mfem::Vector inner(direction);
                    inner *= radius / 2.0;
                    mfem::Vector expected(inner);
                    expected /= std::sqrt(3.0);
                    mfem::Vector mapped = TransformCopy(inner, *cfg, 1);
                    mapped -= expected;
                    EXPECT_LT(mapped.Norml2(), 2.0e-14);
                    for (double interfaceRadius : {radius / 2.0, radius}) {
                        mfem::Vector inside(direction), outside(direction);
                        inside *= interfaceRadius * (1.0 - 1.0e-8);
                        outside *= interfaceRadius * (1.0 + 1.0e-8);
                        mapped = TransformCopy(inside, *cfg, 1);
                        mapped -= TransformCopy(outside, *cfg, interfaceRadius == radius ? 2 : 1);
                        EXPECT_LT(mapped.Norml2(), 1.0e-7 * radius);
                    }
                    mfem::Vector coreInterface(direction);
                    coreInterface *= radius;
                    EXPECT_NEAR(TransformCopy(coreInterface, *cfg, 1).Norml2(), radius, 2.0e-14);
                }
            }
        }
    }
    auto mesh = stroid::GenerateMesh(*cfg);
    ASSERT_NE(mesh.mesh, nullptr);
    ExpectCoreFaceContinuity(*mesh.mesh, 1);
}

TEST_F(stroidTest, MultiBlockCore_ClosedGridSignedJacobiansAndSvdAcrossOrdersAndRefinements) {
    for (int order = 1; order <= 6; ++order) {
        for (int refinement = 0; refinement <= 2; ++refinement) {
            SCOPED_TRACE("order=" + std::to_string(order) + " refinement=" + std::to_string(refinement));
            auto cfg = MultiBlockConfiguration(order, refinement, false);
            auto mesh = stroid::GenerateMesh(*cfg);
            ASSERT_NE(mesh.mesh, nullptr);
            const int factor = 1 << (3 * refinement);
            EXPECT_EQ(mesh.mesh->GetNE(), 13 * factor);
            EXPECT_EQ(CountVolumeAttributes(*mesh.mesh).at(1), 7 * factor);
            ExpectClosedGridCoreConditioning(*mesh.mesh, 1);
        }
    }
}

TEST_F(stroidTest, MultiBlockCore_MapIsScaleInvariantBelowLegacyRadiusCutoff) {
    constexpr double scale = 1.0e-15;
    auto reference = MultiBlockConfiguration(2, 0, true);
    reference->mutate([](stroid::config::MeshConfig& value) { value.r_infinity = 5.0; });
    auto scaled = MultiBlockConfiguration(2, 0, true);
    scaled->mutate([](stroid::config::MeshConfig& value) {
        value.r_core = 2.5e-16;
        value.r_star = 1.0e-15;
        value.r_infinity = 5.0e-15;
    });
    const std::array<std::array<double, 3>, 10> points{{
        {{0.0, 0.0, 0.0}},
        {{0.05, -0.04, 0.1}},
        {{0.125, 0.08, -0.02}},
        {{0.18, -0.09, 0.12}},
        {{-0.2, -0.2, -0.2}},
        {{0.25, 0.12, -0.2}},
        {{0.6, -0.2, 0.4}},
        {{1.0, 0.7, -0.3}},
        {{3.0, -1.3, 0.4}},
        {{-5.0, 2.1, -1.0}}
    }};
    for (const auto& coordinates : points) {
        mfem::Vector point(3);
        for (int component = 0; component < 3; ++component) point(component) = coordinates[component];
        const double logicalRadius = std::max({std::abs(point(0)), std::abs(point(1)), std::abs(point(2))});
        const int attribute = logicalRadius <= 0.25 ? 1 : logicalRadius <= 1.0 ? 2 : 3;
        const auto expected = TransformCopy(point, *reference, attribute);
        point *= scale;
        auto actual = TransformCopy(point, *scaled, attribute);
        actual /= scale;
        for (int component = 0; component < 3; ++component) {
            EXPECT_NEAR(actual(component), expected(component), 2.0e-13)
                << "logical radius=" << logicalRadius << " component=" << component;
        }
    }
}

TEST_F(stroidTest, MultiBlockCore_FlatteningCustomIdsAndExteriorCoordinateRemainConsistent) {
    auto cfg = MultiBlockConfiguration(3, 1, true, 0.2);
    cfg->mutate([](stroid::config::MeshConfig& value) {
        value.core_id = 11;
        value.envelope_id = 17;
        value.vacuum_id = 23;
        value.surface_bdr_id = 31;
        value.inf_bdr_id = 37;
    });
    auto mesh = stroid::GenerateMesh(*cfg);
    ASSERT_NE(mesh.mesh, nullptr);
    const auto volume = CountVolumeAttributes(*mesh.mesh);
    EXPECT_EQ(volume.at(11), 7 * 8);
    EXPECT_EQ(volume.at(17), 6 * 8);
    EXPECT_EQ(volume.at(23), 6 * 8);
    const auto boundary = CountBoundaryAttributes(*mesh.mesh);
    EXPECT_EQ(boundary.at(31), 6 * 4);
    EXPECT_EQ(boundary.at(37), 6 * 4);
    ExpectClosedGridCoreConditioning(*mesh.mesh, 11);
    ExpectCoreFaceContinuity(*mesh.mesh, 11);
    ExpectExteriorCoordinateRange(mesh);
    ExpectExteriorCoordinateBoundaryTraces(mesh);
    mfem::Vector point(3);
    point(0) = 0.25;
    point(1) = 0.25;
    point(2) = 0.25;
    auto mapped = TransformCopy(point, *cfg, 11);
    mapped(2) /= 0.8;
    EXPECT_NEAR(mapped.Norml2(), 0.25, 2.0e-14);
}

TEST_F(stroidTest, MultiBlockCore_OuterMappingAndSignedStellarVolumeMatchLegacy) {
    auto cfg = MultiBlockConfiguration(3, 1, true);
    auto legacyCfg = MultiBlockConfiguration(3, 1, true);
    legacyCfg->mutate([](stroid::config::MeshConfig& value) { value.core_mapping = "spherified"; });
    const double coreRadius = (*cfg)->r_core.value();
    const double stellarRadius = (*cfg)->r_star.value();
    const double infinityRadius = (*cfg)->r_infinity.value();
    for (int axis = 0; axis < 3; ++axis) {
        for (double sign : {-1.0, 1.0}) {
            for (double a : {-1.0, -0.3, 0.0, 0.8, 1.0}) {
                for (double b : {-1.0, 0.0, 0.4, 1.0}) {
                    mfem::Vector direction(3);
                    direction(axis) = sign;
                    direction((axis + 1) % 3) = a;
                    direction((axis + 2) % 3) = b;
                    for (double radius : {coreRadius, (coreRadius + stellarRadius) / 2.0, stellarRadius,
                                          (stellarRadius + infinityRadius) / 2.0, infinityRadius}) {
                        mfem::Vector point(direction);
                        point *= radius;
                        const int attribute = radius <= stellarRadius ? 2 : 3;
                        auto difference = TransformCopy(point, *cfg, attribute);
                        difference -= TransformCopy(point, *legacyCfg, attribute);
                        EXPECT_LT(difference.Norml2(), 2.0e-14 * infinityRadius);
                    }
                }
            }
        }
    }
    auto mesh = stroid::GenerateMesh(*cfg);
    auto legacy = stroid::GenerateMesh(*legacyCfg);
    const auto signedStellarVolume = [](mfem::Mesh& candidate) {
        double volume = 0.0;
        for (int element = 0; element < candidate.GetNE(); ++element) {
            if (candidate.GetAttribute(element) == 3) continue;
            auto* transformation = candidate.GetElementTransformation(element);
            const auto& rule = mfem::IntRules.Get(transformation->GetGeometryType(), 3 * transformation->Order() + 2);
            for (int q = 0; q < rule.GetNPoints(); ++q) {
                const auto& point = rule.IntPoint(q);
                transformation->SetIntPoint(&point);
                volume += point.weight * transformation->Jacobian().Det();
            }
        }
        return volume;
    };
    const double newVolume = signedStellarVolume(*mesh.mesh);
    const double oldVolume = signedStellarVolume(*legacy.mesh);
    EXPECT_GT(newVolume, 0.0);
    EXPECT_NEAR(newVolume, oldVolume, 2.0e-11 * oldVolume);
}

TEST_F(stroidTest, MultiBlockCore_SaveLoadConfigAndRefinementPreserveContracts) {
    for (const bool external : {false, true}) {
        SCOPED_TRACE(external);
        auto cfg = MultiBlockConfiguration(3, 0, external);
        auto original = stroid::GenerateMesh(*cfg);
        EXPECT_EQ(original.type, stroid::MFEM_MESH_TYPE::SERIAL);
        const auto path = std::filesystem::temp_directory_path() /
            (external ? "stroid_multiblock_external_round_trip.smesh" : "stroid_multiblock_stellar_round_trip.smesh");
        stroid::IO::SaveStroidMesh(original, path.string(), "Multi-block core regression");
        auto result = stroid::IO::LoadStroidMesh(path.string());
        ASSERT_TRUE(result.has_value()) << result.error();
        auto loaded = std::move(*result);
        EXPECT_EQ(loaded.type, stroid::MFEM_MESH_TYPE::SERIAL);
        ASSERT_NE(loaded.mesh, nullptr);
        ASSERT_NE(loaded.reference_mesh, nullptr);
        EXPECT_EQ(loaded.config.core_mapping.value(), "multi_block");
        EXPECT_EQ(loaded.config.include_external_domain.value(), external);
        EXPECT_EQ(loaded.mesh->GetNE(), original.mesh->GetNE());
        ASSERT_EQ(loaded.mesh->GetNodes()->Size(), original.mesh->GetNodes()->Size());
        for (int dof = 0; dof < original.mesh->GetNodes()->Size(); ++dof) {
            EXPECT_NEAR((*loaded.mesh->GetNodes())(dof), (*original.mesh->GetNodes())(dof), 2.0e-14);
        }
        stroid::refinement::UniformRefinement(loaded, 1);
        EXPECT_EQ(loaded.refinement_levels, 1);
        EXPECT_EQ(loaded.mesh->GetNE(), original.mesh->GetNE() * 8);
        EXPECT_EQ(CountVolumeAttributes(*loaded.mesh).at(1), 7 * 8);
        ExpectClosedGridCoreConditioning(*loaded.mesh, 1);
        ExpectCoreFaceContinuity(*loaded.mesh, 1);
        if (external) {
            ExpectExteriorCoordinateRange(loaded);
            ExpectExteriorCoordinateBoundaryTraces(loaded);
        } else {
            EXPECT_EQ(loaded.exterior_coordinate, nullptr);
        }
        std::error_code error;
        std::filesystem::remove(path, error);
        EXPECT_FALSE(error);
    }

    auto legacyCfg = MultiBlockConfiguration(2, 0, false);
    legacyCfg->mutate([](stroid::config::MeshConfig& value) { value.core_mapping = "spherified"; });
    auto legacy = stroid::GenerateMesh(*legacyCfg);
    EXPECT_EQ(legacy.type, stroid::MFEM_MESH_TYPE::SERIAL);
    const auto path = std::filesystem::temp_directory_path() / "stroid_core_mapping_legacy_round_trip.smesh";
    stroid::IO::SaveStroidMesh(legacy, path.string(), "Legacy core mapping default regression");
    std::ifstream input(path);
    std::string contents((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    const auto marker = contents.find("\ncore_mapping:");
    ASSERT_NE(marker, std::string::npos);
    const auto fieldStart = marker + 1;
    const auto newline = contents.find('\n', fieldStart);
    ASSERT_NE(newline, std::string::npos);
    contents.erase(fieldStart, newline - fieldStart + 1);
    std::istringstream legacyStream(contents);
    auto restored = stroid::IO::ParseStroidMesh(legacyStream);
    ASSERT_TRUE(restored.has_value()) << restored.error();
    EXPECT_EQ(restored->config.core_mapping.value(), "spherified");
    EXPECT_EQ(CountVolumeAttributes(*restored->mesh).at(1), 1);
    std::error_code error;
    std::filesystem::remove(path, error);
    EXPECT_FALSE(error);
}
