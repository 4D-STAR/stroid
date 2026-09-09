#include "stroid/topology/curvilinear.h"
#include "stroid/topology/mapping.h"

#include <iostream>

namespace {
    double compute_exterior_coordinate(
        const mfem::Vector& logical_position,
        const int attribute,
        const fourdst::config::Config<stroid::config::MeshConfig>& config
    ) {
        if (!config->include_external_domain.value_or(true) || attribute != static_cast<int>(config->vacuum_id.value_or(3))) return 0.0;

        const double r_star = config->r_star.value_or(1.0);
        const double r_infinity = config->r_infinity.value_or(6.0);
        const double radial_extent = r_infinity - r_star;

        if (!std::isfinite(r_star) || !std::isfinite(r_infinity) || r_star <= 0.0 || radial_extent <= 0.0) {
            throw std::invalid_argument("Exterior-coordinate construction requires finite radii with 0 < r_star < r_infinity.");
        }

        double logical_radius = 0.0;
        for (int d = 0; d < logical_position.Size(); ++d) {
            if (!std::isfinite(logical_position(d))) throw std::runtime_error("Reference mesh produced a non-finite logical position.");
            logical_radius = std::max(logical_radius, std::abs(logical_position(d)));
        }

        double coordinate = (logical_radius - r_star) / radial_extent;
        const double tolerance = std::min(1.0e-8, 1024.0 * std::numeric_limits<double>::epsilon() *
            std::max(std::abs(r_star), std::abs(r_infinity)) / radial_extent);

        if (coordinate < -tolerance || coordinate > 1.0 + tolerance) {
            throw std::runtime_error(std::format("Logical exterior coordinate {} lies outside [0, 1].", coordinate));
        }

        if (std::abs(coordinate) <= tolerance) coordinate = 0.0;
        if (std::abs(coordinate - 1.0) <= tolerance) coordinate = 1.0;
        return coordinate;
    }
}

namespace stroid::topology {
    void PromoteToHighOrder(mfem::Mesh &mesh, const fourdst::config::Config<config::MeshConfig> &config) {
        mesh.SetCurvature(config->order.value(), false, mesh.SpaceDimension(), mfem::Ordering::byNODES);
    }

    void ProjectMesh(mfem::Mesh &mesh, const fourdst::config::Config<config::MeshConfig> &config) {
        if (!mesh.GetNodes()) {
            std::cerr << "Error: Mesh has no nodes to project. Call PromoteToHighOrder first." << std::endl;
            return;
        }

        mfem::GridFunction& nodes = *mesh.GetNodes(); // Already confirmed not null
        const mfem::FiniteElementSpace* fes = nodes.FESpace();

        const int vDim = fes->GetVDim();
        const int nDofs = fes->GetNDofs();
        const int nElem = mesh.GetNE();

        std::vector<bool> processed(nDofs, false);
        mfem::Array<int> dofs;
        mfem::Vector pos(vDim);

        for (int elemID = 0; elemID < nElem; ++elemID) {
            const int attrID = mesh.GetAttribute(elemID);
            fes->GetElementDofs(elemID, dofs);

            for (int dofID = 0; dofID < dofs.Size(); ++dofID) {
                const int scalar_dof = dofs[dofID] >= 0 ? dofs[dofID] : -1 - dofs[dofID];

                if (processed[scalar_dof]) {
                    continue; // Skip already processed dofs. This avoids doing multiple transformations of a node if it was already transformed by a neighbor
                }

                for (int d = 0; d < vDim; ++d) {
                    pos(d) = nodes(fes->DofToVDof(scalar_dof, d));
                }

                TransformPoint(pos, config, attrID);

                for (int d = 0; d < vDim; ++d) {
                    nodes(fes->DofToVDof(scalar_dof, d)) = pos(d);
                }

                processed[scalar_dof] = true;
            }
        }

        // A mapped hanging node must lie on the coarse polynomial face. Mapping
        // every node independently does not preserve this geometric constraint.
        mfem::Vector true_nodes;
        nodes.GetTrueDofs(true_nodes);
        nodes.SetFromTrueDofs(true_nodes);
        mesh.NodesUpdated();
    }

    std::unique_ptr<ScalarMeshField> BuildExteriorCoordinate(
        mfem::Mesh& mesh,
        mfem::Mesh& reference_mesh,
        const fourdst::config::Config<config::MeshConfig>& config
    ) {
        if (!config->include_external_domain.value_or(true)) return nullptr;
        if (mesh.Dimension() != reference_mesh.Dimension() || mesh.SpaceDimension() != reference_mesh.SpaceDimension()) {
            throw std::invalid_argument("Primary and reference meshes must have matching dimensions when constructing the exterior coordinate.");
        }
        if (mesh.GetNE() != reference_mesh.GetNE()) {
            throw std::invalid_argument("Primary and reference meshes must have the same number of elements when constructing the exterior coordinate.");
        }
        if (mesh.GetNodalFESpace() == nullptr) {
            throw std::invalid_argument("Exterior-coordinate construction requires a primary mesh with a nodal finite-element space.");
        }

        for (int element_id = 0; element_id < mesh.GetNE(); ++element_id) {
            if (mesh.GetElementGeometry(element_id) != reference_mesh.GetElementGeometry(element_id)) {
                throw std::invalid_argument(std::format("Primary and reference element {} have different geometries.", element_id));
            }
            if (mesh.GetAttribute(element_id) != reference_mesh.GetAttribute(element_id)) {
                throw std::invalid_argument(std::format("Primary and reference element {} have different attributes.", element_id));
            }
        }

        auto field = std::make_unique<ScalarMeshField>();
        const mfem::FiniteElementCollection* collection = mesh.GetNodalFESpace()->FEColl();
        field->space = std::make_unique<mfem::FiniteElementSpace>(&mesh, collection);
        field->values = std::make_unique<mfem::GridFunction>(field->space.get());
        *field->values = 0.0;

        const int scalar_dofs = field->space->GetNDofs();
        std::vector<bool> processed(static_cast<size_t>(scalar_dofs), false);
        mfem::Array<int> element_dofs;
        mfem::Vector logical_position(reference_mesh.SpaceDimension());

        const double consistency_tolerance = 4096.0 * std::numeric_limits<double>::epsilon();

        for (int element_id = 0; element_id < mesh.GetNE(); ++element_id) {
            const mfem::FiniteElement& element = *field->space->GetFE(element_id);
            const mfem::IntegrationRule& nodes = element.GetNodes();
            mfem::ElementTransformation* reference_transformation = reference_mesh.GetElementTransformation(element_id);

            if (reference_transformation == nullptr) throw std::runtime_error(std::format("Reference element {} has no element transformation.", element_id));

            field->space->GetElementDofs(element_id, element_dofs);
            if (nodes.GetNPoints() != element_dofs.Size()) {
                throw std::runtime_error(std::format("Element {} has {} nodal points but {} scalar DOFs.", element_id, nodes.GetNPoints(), element_dofs.Size()));
            }

            for (int local_dof = 0; local_dof < element_dofs.Size(); ++local_dof) {
                const int encoded_dof = element_dofs[local_dof];
                const int global_dof = encoded_dof >= 0 ? encoded_dof : -1 - encoded_dof;

                if (global_dof < 0 || global_dof >= scalar_dofs) {
                    throw std::runtime_error(std::format("Element {} references invalid scalar DOF {}.", element_id, global_dof));
                }

                reference_transformation->Transform(nodes.IntPoint(local_dof), logical_position);
                const double coordinate = compute_exterior_coordinate(logical_position, mesh.GetAttribute(element_id), config);

                if (processed[static_cast<size_t>(global_dof)]) {
                    const double existing_coordinate = (*field->values)(global_dof);
                    if (std::abs(existing_coordinate - coordinate) > consistency_tolerance) {
                        throw std::runtime_error(std::format("Exterior coordinate is inconsistent at shared scalar DOF {}: existing value {}, new value {} from element {}.", global_dof, existing_coordinate, coordinate, element_id));
                    }
                    continue;
                }

                (*field->values)(global_dof) = coordinate;
                processed[static_cast<size_t>(global_dof)] = true;
            }
        }

        mfem::Vector true_values;
        field->values->GetTrueDofs(true_values);
        field->values->SetFromTrueDofs(true_values);

        for (int dof = 0; dof < scalar_dofs; ++dof) {
            if (!processed[static_cast<size_t>(dof)]) throw std::runtime_error(std::format("Exterior-coordinate scalar DOF {} was not assigned.", dof));
            const double coordinate = (*field->values)(dof);
            if (!std::isfinite(coordinate) || coordinate < 0.0 || coordinate > 1.0) {
                throw std::runtime_error(std::format("Exterior-coordinate scalar DOF {} has invalid value {}.", dof, coordinate));
            }
        }

        return field;
    }
}
