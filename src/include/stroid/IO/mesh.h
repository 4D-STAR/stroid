#pragma once
#include <string>
#include "mfem.hpp"

namespace stroid::IO {
    /**
     * @brief Visualization modes for GLVis display.
     */
    enum class VISUALIZATION_MODE : uint8_t {
        /** @brief No attribute visualization (default rendering). */
        NONE,
        /** @brief Color elements by their element attribute/ID. */
        ELEMENT_ID,
        /** @brief Color boundary-adjacent elements by boundary attribute/ID. */
        BOUNDARY_ELEMENT_ID
    };

    /**
     * @brief Save a mesh to MFEM's native `.mesh` format.
     * @param mesh Mesh to serialize.
     * @param filename Output path (including extension).
     */
    void SaveMesh(const mfem::Mesh& mesh, const std::string& filename);
    /**
     * @brief Save a mesh as a ParaView VTU dataset.
     * @param mesh Mesh to export.
     * @param exportName Output base name (ParaView will add extensions).
     */
    void SaveVTU(mfem::Mesh& mesh, const std::string& exportName);
    /**
     * @brief Stream a mesh to a running GLVis server for interactive viewing.
     * @param mesh Mesh to display.
     * @param title Window title shown in GLVis.
     * @param mode Attribute visualization mode.
     * @param vishost GLVis server host.
     * @param visport GLVis server port.
     */
    void ViewMesh(mfem::Mesh &mesh, const std::string& title, VISUALIZATION_MODE mode, const std::string &vishost, int visport);
    /**
     * @brief Visualize boundary face valence (1=surface, 2=internal).
     * @param mesh Mesh whose boundary faces are inspected.
     */
    void VisualizeFaceValence(mfem::Mesh& mesh);
}