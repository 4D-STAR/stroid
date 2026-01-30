#include <gtest/gtest.h>

#include "fourdst/config/config.h"
#include "stroid/config/config.h"
#include "stroid/IO/mesh.h"
#include "stroid/topology/curvilinear.h"
#include "stroid/topology/mapping.h"
#include "stroid/topology/topology.h"

#include <cmath>
#include <filesystem>
#include <cstdlib>
#include <string>

namespace {

constexpr double kPi = 3.14159265358979323846;

using Config = fourdst::config::Config<stroid::config::MeshConfig>;


std::filesystem::path GetSourceRoot() {
    if (const char* env = std::getenv("MESON_SOURCE_ROOT")) {
        return std::filesystem::path(env);
    }
    return std::filesystem::current_path();
}

Config LoadConfigFromRepo(const std::filesystem::path& relative_path) {
    Config cfg;
    cfg.load((GetSourceRoot() / relative_path).string());
    return cfg;
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

} // namespace

/**
 * @brief Test suite for the Stroid library
 */
class stroidTest : public ::testing::Test {};

TEST_F(stroidTest, BuildSkeleton_DefaultCounts) {
    const Config cfg;
    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);

    ASSERT_NE(mesh, nullptr);
    EXPECT_EQ(mesh->Dimension(), 3);
    EXPECT_EQ(mesh->GetNV(), 16);
    EXPECT_EQ(mesh->GetNE(), 7);
    EXPECT_EQ(mesh->GetNBE(), 6);
}


TEST_F(stroidTest, Finalize_RefinementIncreasesElements) {
    const Config cfg;
    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    const int initial_elements = mesh->GetNE();

    stroid::topology::Finalize(*mesh, cfg);

    EXPECT_GT(mesh->GetNE(), initial_elements);
    EXPECT_TRUE(mesh->Conforming());
}

TEST_F(stroidTest, PromoteToHighOrder_SetsNodes) {
    const Config cfg;
    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    stroid::topology::Finalize(*mesh, cfg);

    stroid::topology::PromoteToHighOrder(*mesh, cfg);

    EXPECT_NE(mesh->GetNodes(), nullptr);
    EXPECT_TRUE(IsFiniteMeshNodes(*mesh));
}

TEST_F(stroidTest, ProjectMesh_ProducesFiniteNodes) {
    const Config cfg;
    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    stroid::topology::Finalize(*mesh, cfg);
    stroid::topology::PromoteToHighOrder(*mesh, cfg);

    stroid::topology::ProjectMesh(*mesh, cfg);

    EXPECT_TRUE(IsFiniteMeshNodes(*mesh));
}

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

TEST_F(stroidTest, ApplySpheroidal_FlattensZ) {
    const Config cfg = LoadConfigFromRepo("configs/test_flattening.toml");

    mfem::Vector pos(3);
    pos(0) = 0.0;
    pos(1) = 0.0;
    pos(2) = 10.0;

    stroid::topology::ApplySpheroidal(pos, cfg);

    EXPECT_NEAR(pos(2), 8.0, 1e-12);
}

TEST_F(stroidTest, ApplyKelvin_ExpandsOutsideStar) {
    const Config cfg;

    mfem::Vector pos(3);
    pos(0) = 5.5;
    pos(1) = 0.0;
    pos(2) = 0.0;

    stroid::topology::ApplyKelvin(pos, cfg);

    EXPECT_NEAR(pos(0), 6.0, 1e-12);
    EXPECT_NEAR(pos(1), 0.0, 1e-12);
    EXPECT_NEAR(pos(2), 0.0, 1e-12);
}

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

TEST_F(stroidTest, TransformPoint_AxisOutsideStar_KelvinExpands) {
    const Config cfg;

    mfem::Vector pos(3);
    pos(0) = 5.5;
    pos(1) = 0.0;
    pos(2) = 0.0;

    stroid::topology::TransformPoint(pos, cfg, 0);

    EXPECT_NEAR(pos(0), 6.0, 1e-12);
    EXPECT_NEAR(pos(1), 0.0, 1e-12);
    EXPECT_NEAR(pos(2), 0.0, 1e-12);
}

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
