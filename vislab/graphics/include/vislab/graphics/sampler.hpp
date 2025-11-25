#pragma once

#include <vislab/core/object.hpp>

#include <Eigen/Eigen>

namespace vislab
{
    /**
     * @brief Base class for random number generators.
     */
    class Sampler : public Object
    {
        VISLAB_INTERFACE(Sampler, Object)

    public:
        /**
         * @brief Seeds the random number sequence.
         * @param seed_offset Seed offset of the sequence.
         */
        virtual void seed(uint64_t seed_offset);

        /**
         * @brief Advances the sample index.
         */
        void advance();

        /**
         * @brief Gets the sample index.
         * @return Sample index.
         */
        uint32_t getSampleIndex() const;

        /**
         * @brief Gets the total number of samples.
         * @return Total number of samples.
         */
        uint32_t getSampleCount() const;

        /**
         * @brief Gets the base seed of the sequence.
         * @return Base seed of sequence.
         */
        uint64_t getBaseSeed() const;

        /**
         * @brief Sets the base seed of the sequence.
         * @param baseSeed Base seed of sequence.
         */
        void setBaseSeed(uint64_t baseSeed);

        /**
         * @brief Generates a random number between [0,1].
         * @return Random number.
         */
        virtual double next_1d() = 0;

        /**
         * @brief Generates two random numbers between [0,1].
         * @return Two random numbers.
         */
        virtual Eigen::Vector2d next_2d() = 0;

    protected:
        /**
         * @brief Constructor.
         */
        Sampler();

        /**
         * @brief Copy-Constructor.
         * @param other Sampler to copy from.
         */
        Sampler(const Sampler& other);

        /**
         * @brief Base seed of the sequence.
         */
        uint64_t mBaseSeed;

        /**
         * @brief Number of samples in the sequence.
         */
        uint32_t mSampleCount;

        /**
         * @brief Index of the current sample.
         */
        uint32_t mSampleIndex;
    };
}
