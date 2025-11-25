#pragma once

#include <stdint.h>

namespace vislab
{
    /**
     * @brief Generate pseudo random numbers using the Tiny Encryption Algorithm (TEA) by David Wheeler and Roger Needham. See http://dx.doi.org/10.2312/EGGH/HPG10/133-141
     */
    class TinyEncryptionAlgorithm
    {
        /**
         * @brief Generate a 32-bit unsigned int.
         * @param v0
         * @param v1
         * @param rounds
         * @return
         */
        static uint32_t sample_tea_32(uint32_t v0, uint32_t v1, int rounds = 4);

        /**
         * @brief Generate a 64-bit unsigned int.
         * @param v0
         * @param v1
         * @param rounds
         * @return
         */
        static uint64_t sample_tea_64(uint32_t v0, uint32_t v1, int rounds = 4);

        /**
         * @brief Generate a 32-bit float.
         * @param v0
         * @param v1
         * @param rounds
         * @return
         */
        static float sample_tea_float32(uint32_t v0, uint32_t v1, int rounds = 4);

        /**
         * @brief Generate a 64-bit double.
         * @param v0
         * @param v1
         * @param rounds
         * @return
         */
        static double sample_tea_float64(uint32_t v0, uint32_t v1, int rounds = 4);
    };
}
