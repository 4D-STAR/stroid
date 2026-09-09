#pragma once
#include <string>
#include <expected>
#include <istream>
#include <memory>

#include "mfem.hpp"

#include "stroid/utils/types.h"

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

    void SaveStroidMesh(const StroidMesh& mesh, const std::string& filename, const std::string& comment="");

    /**
     * @brief Save a mesh to MFEM's native `.mesh` format.
     * @param mesh Mesh to serialize.
     * @param filename Output path (including extension).
     */
    void SaveMesh(const mfem::Mesh& mesh, const std::string& filename);

    /**
     * @brief Overload of SaveMesh which accepts a StroidMesh type and will internally unpack it
     * @param mesh StroidMesh to serialize.
     * @param filename Path to save to
     *
     * @note This function is a utility wrapper to save a StroidMesh object in MFEM's native .mesh format. Data other than the mesh pointer
     * in StroidMesh **will not be saved** (e.g. the reference mesh, the number of refinement levels, etc..). If you need to serialize an
     * entire StroidMesh then please use the stroid::IO::SaveStroidMesh function
     */
    void SaveMesh(const stroid::StroidMesh& mesh, const std::string& filename);
    /**
     * @brief Save a mesh as a ParaView VTU dataset.
     * @param mesh Mesh to export.
     * @param exportName Output base name (ParaView will add extensions).
     */
    void SaveVTU(mfem::Mesh& mesh, const std::string& exportName);

    /**
     * @brief Overload of SaveVTU which accepts a StroidMesh type and will internally unpack it
     * @param mesh StroidMesh to serialize.
     * @param filename Path to save to
     *
     * @note This function is a utility wrapper to save a StroidMesh object in MFEM's native .mesh format. Data other than the mesh pointer
     * in StroidMesh **will not be saved** (e.g. the reference mesh, the number of refinement levels, etc..). If you need to serialize an
     * entire StroidMesh then please use the stroid::IO::SaveStroidVTU function
     */
    void SaveVTU(const stroid::StroidMesh& mesh, const std::string& exportName);

    /**
     * @brief Make a display-only mesh copy with matching face subdivisions.
     * This is purely for visualization and should not be used for any science goals.
     * @param mesh Source mesh whose geometry is to be displayed.
     * @return Independently owned copy with no hanging faces.
     */
    std::unique_ptr<mfem::Mesh> MakeConformingVisualizationMesh(const mfem::Mesh& mesh);

    /**
     * @brief Stream a mesh to a running GLVis server for interactive viewing.
     * @param mesh Mesh to display.
     * @param title Window title shown in GLVis.
     * @param mode Attribute visualization mode.
     * @param vishost GLVis server host.
     * @param visport GLVis server port.
     * @param conforming_display Refine a display-only copy at hanging interfaces
     * to prevent GLVis tessellation gaps. Set false to inspect the original
     * element layout, which can show rendering gaps on curved interfaces.
     *
     */
    void ViewMesh(mfem::Mesh &mesh, const std::string& title, VISUALIZATION_MODE mode, const std::string &vishost, int visport, bool conforming_display=true);

    void ViewMesh(const stroid::StroidMesh& mesh, const std::string& title, VISUALIZATION_MODE mode, const std::string &vishost, int visport, bool conforming_display=true);

    /**
     * @brief Color boundary-adjacent elements by face valence (1=surface, 2=internal).
     * Untagged elements are zero; elements touching several tagged faces use
     * the maximum valence. Values are computed before display subdivision.
     * @param mesh Mesh whose boundary faces are inspected.
     * @param vishost GLVis server host.
     * @param visport GLVis server port.
     * @param conforming_display Use the same display-only refinement as ViewMesh.
     */
    void VisualizeFaceValence(mfem::Mesh& mesh, const std::string &vishost, int visport, bool conforming_display=true);

    void VisualizeFaceValence(const stroid::StroidMesh& mesh, const std::string &vishost, int visport, bool conforming_display=true);

    std::expected<StroidMesh, std::string> ParseStroidMesh(std::istream& is);
    std::expected<StroidMesh, std::string> LoadStroidMesh(const std::string& filename);

#ifdef MFEM_USE_MPI
    std::expected<StroidMesh, std::string> ParseStroidMesh(std::istream& is, MPI_Comm comm);
    std::expected<StroidMesh, std::string> LoadStroidMesh(const std::string& filename, MPI_Comm comm);
#endif
}
