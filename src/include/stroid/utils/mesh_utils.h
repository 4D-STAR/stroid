#pragma once

#include "mfem.hpp"

namespace stroid::utils {
    /**
     * @brief Mark elements with negative Jacobian determinant.
     * @details Elements detected as flipped are assigned attribute 999.
     * @param mesh Mesh to scan and update in-place.
     */
    void MarkFlippedElements(mfem::Mesh& mesh);
    /**
     * @brief Mark boundary elements whose outward normal orientation is flipped.
     * @details Boundary elements detected as flipped are assigned attribute 500.
     * @param mesh Mesh to scan and update in-place.
     */
    void MarkFlippedBoundaryElements(mfem::Mesh& mesh);
}