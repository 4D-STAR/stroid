#include "mfem.hpp"

#include "stroid/refinement/uniform.h"
#include "stroid/utils/types.h"
#include "stroid/utils/mesh_utils.h"
#include "stroid/exceptions/exceptions.h"
#include "stroid/topology/curvilinear.h"
#include "stroid/topology/topology.h"
#include "stroid/topology/optimize.h"

#include <limits>

namespace stroid::refinement {
    void UniformRefinement(StroidMesh &mesh, const size_t levels) {
        if (!mesh.reference_mesh) {
            throw exceptions::StroidMissingReferenceMesh("UniformRefinement requires a reference mesh to be present in the StroidMesh object. This should be present by construction and the fact that is is missing represents a bug. Please report this to the stroid developers on GitHub or by email at emily.boudreaux@dartmouth.edu");
        }

        if (levels == 0) {
            return;
        }

        if (!mesh.mesh) {
            throw exceptions::StroidMissingReferenceMesh(
                "UniformRefinement requires a primary mesh to be present in the StroidMesh object. This should be present by construction and the fact that it is missing represents a bug. Please report this to the stroid developers on GitHub or by email at emily.boudreaux@dartmouth.edu");
        }
        if (levels > std::numeric_limits<size_t>::max() - mesh.refinement_levels) {
            throw std::overflow_error("Uniform refinement level count would overflow.");
        }

        StroidMesh refined;
        refined.type = mesh.type;
        refined.config = mesh.config;
        refined.refinement_levels = mesh.refinement_levels + levels;
        refined.reference_mesh = std::make_unique<mfem::Mesh>(*mesh.reference_mesh);
        for (size_t i = 0; i < levels; i++) {
            refined.reference_mesh->UniformRefinement();
        }

        fourdst::config::Config<config::MeshConfig> cfg;
        auto Mutator = [&mesh](config::MeshConfig& orig) {
            orig = mesh.config;
        };

        cfg.mutate(Mutator);

        refined.mesh = utils::BuildProjected(*refined.reference_mesh, cfg);
        topology::OptimizeMesh(*refined.mesh, cfg);
        refined.exterior_coordinate = topology::BuildExteriorCoordinate(*refined.mesh, *refined.reference_mesh, cfg);

        mesh.mesh.swap(refined.mesh);
        mesh.reference_mesh.swap(refined.reference_mesh);
        mesh.exterior_coordinate.swap(refined.exterior_coordinate);
        mesh.refinement_levels = refined.refinement_levels;
    }
}
