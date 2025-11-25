#pragma once

#include <stdint.h>

namespace vislab
{
    /**
     * \brief Flags used to classify the different types of lobes that are implemented in a BSDF instance.
     */
    enum class EBSDFFlag : uint32_t
    {
        // =============================================================
        //                      BSDF lobe types
        // =============================================================

        /**
         * @brief No flags set (default value)
         */
        Empty = 0x00000,

        /**
         * @brief 'null' scattering event, i.e. particles do not undergo deflection
         */
        Null = 0x00001,

        /**
         * @brief Ideally diffuse reflection
         */
        DiffuseReflection = 0x00002,

        /**
         * @brief Ideally diffuse transmission
         */
        DiffuseTransmission = 0x00004,

        /**
         * @brief Glossy reflection
         */
        GlossyReflection = 0x00008,

        /**
         * @brief Glossy transmission
         */
        GlossyTransmission = 0x00010,

        /**
         * @brief Reflection into a discrete set of directions
         */
        DeltaReflection = 0x00020,

        /**
         * @brief Transmission into a discrete set of directions
         */
        DeltaTransmission = 0x00040,

        /**
         * @brief Reflection into a 1D space of directions
         */
        Delta1DReflection = 0x00080,

        /**
         * @brief Transmission into a 1D space of directions
         */
        Delta1DTransmission = 0x00100,

        // =============================================================
        //!                   Other lobe attributes
        // =============================================================

        /**
         * @brief The lobe is not invariant to rotation around the normal
         */
        Anisotropic = 0x01000,

        /**
         * @brief The BSDF depends on the UV coordinates
         */
        SpatiallyVarying = 0x02000,

        /**
         * @brief Flags non-symmetry (e.g. transmission in dielectric materials)
         */
        NonSymmetric = 0x04000,

        /**
         * @brief Supports interactions on the front-facing side
         */
        FrontSide = 0x08000,

        /**
         * @brief Supports interactions on the back-facing side
         */
        BackSide = 0x10000,

        /**
         * @brief Does the implementation require access to texture-space differentials
         */
        NeedsDifferentials = 0x20000,

        // =============================================================
        //!                 Compound lobe attributes
        // =============================================================

        /**
         * @brief Any reflection component (scattering into discrete, 1D, or 2D set of directions)
         */
        Reflection = DiffuseReflection | DeltaReflection |
                     Delta1DReflection | GlossyReflection,

        /**
         * @brief Any transmission component (scattering into discrete, 1D, or 2D set of directions)
         */
        Transmission = DiffuseTransmission | DeltaTransmission |
                       Delta1DTransmission | GlossyTransmission | Null,

        /**
         * @brief Diffuse scattering into a 2D set of directions
         */
        Diffuse = DiffuseReflection | DiffuseTransmission,

        /**
         * @brief Non-diffuse scattering into a 2D set of directions
         */
        Glossy = GlossyReflection | GlossyTransmission,

        /**
         * @brief Scattering into a 2D set of directions
         */
        Smooth = Diffuse | Glossy,

        /**
         * @brief Scattering into a discrete set of directions
         */
        Delta = Null | DeltaReflection | DeltaTransmission,

        /**
         * @brief Scattering into a 1D space of directions
         */
        Delta1D = Delta1DReflection | Delta1DTransmission,

        /**
         * @brief Any kind of scattering
         */
        All = Diffuse | Glossy | Delta | Delta1D
    };
}
