#pragma once

#include <stdint.h>

namespace vislab
{
    /**
     * @brief Enumeration of light flags. Renderers might adjust their implementation based on the type.
     */
    enum class ELightFlag : uint32_t
    {
        // =============================================================
        //                      Light types
        // =============================================================

        /**
         * @brief No flags set (default value)
         */
        None = 0x00000,

        /**
         * @brief The light lies at a single point in space
         */
        DeltaPosition = 0x00001,

        /**
         * @brief The light emits light in a single direction
         */
        DeltaDirection = 0x00002,

        /**
         * @brief The light is placed at infinity (e.g. environment maps)
         */
        Infinite = 0x00004,

        /**
         * @brief The light is attached to a surface (e.g. area lights)
         */
        Surface = 0x00008,

        // =============================================================
        //!                   Other lobe attributes
        // =============================================================

        /**
         * @brief The emission depends on the UV coordinates
         */
        SpatiallyVarying = 0x00010,

        // =============================================================
        //!                 Compound lobe attributes
        // =============================================================

        /**
         * @brief Delta function in either position or direction
         */
        Delta = DeltaPosition | DeltaDirection,
    };
}
