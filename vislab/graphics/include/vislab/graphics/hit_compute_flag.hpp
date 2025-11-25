#pragma once

#include <stdint.h>

namespace vislab
{
    /**
     * @brief This list of flags is used to determine which members of SurfaceInteraction should be computed when calling compute_surface_interaction().
     */
    enum class EHitComputeFlag : uint32_t
    {
        // =============================================================
        //             Surface interaction compute flags
        // =============================================================

        /**
         * @brief No flags set
         */
        None = 0x00000,

        /**
         * @brief Compute position and geometric normal
         */
        Minimal = 0x00001,

        /**
         * @brief Compute UV coordinates
         */
        UV = 0x00002,

        /**
         * @brief Compute position partials wrt. UV coordinates
         */
        dPdUV = 0x00004,

        /**
         * @brief Compute shading normal and shading frame
         */
        ShadingFrame = 0x00008,

        // =============================================================
        //                 Compound compute flags
        // =============================================================

        /**
         * @brief Compute all fields of the surface interaction data structure (default)
         */
        All = UV | dPdUV | ShadingFrame,
    };

    /**
     * @brief Tests if f0 and f1 share a common flag.
     * @param f0 First flag collection.
     * @param f1 Second flag collection.
     * @return True if there is at least one common flag.
     */
    constexpr bool hasFlag(EHitComputeFlag f0, EHitComputeFlag f1)
    {
        return ((uint32_t)f0 & (uint32_t)f1) != 0;
    }
}
