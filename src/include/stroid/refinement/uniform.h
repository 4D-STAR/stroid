#pragma once

#include "stroid/utils/types.h"

namespace stroid::refinement {
    /**
     * @brief Refine every current leaf, preserving any existing nonconforming hierarchy.
     * Rebuilds constrained geometry and the exterior coordinate from the refined
     * reference mesh. The saved generation configuration is unchanged; the
     * refinement counter and actual regional depths increase by @p levels.
     */
    void UniformRefinement(StroidMesh& mesh, size_t levels);
}
