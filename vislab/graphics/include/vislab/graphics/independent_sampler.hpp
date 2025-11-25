#pragma once

#include "sampler.hpp"

#include <random>

namespace vislab
{
    /**
     * @brief Independent uniform sampler for random number generation.
     */
    class IndependentSampler : public Sampler
    {
        VISLAB_OBJECT(IndependentSampler, Sampler)

    public:
        /**
         * @brief Constructor.
         */
        IndependentSampler();

        /**
         * @brief Copy-Constructor.
         * @param other Sampler to copy from.
         */
        IndependentSampler(const IndependentSampler& other);

        /**
         * @brief Seeds the random number sequence.
         * @param seed_offset Seed offset.
         */
        void seed(uint64_t seed_offset = std::default_random_engine::default_seed) override;

        /**
         * @brief Generates a random number between [0,1].
         * @return Random number.
         */
        double next_1d() override;

        /**
         * @brief Generates two random numbers between [0,1].
         * @return Two random numbers.
         */
        Eigen::Vector2d next_2d() override;

    private:
        /**
         * @brief Random number generation engine.
         */
        std::default_random_engine mRndGenerator;

        /**
         * @brief Uniform real distribution.
         */
        std::uniform_real_distribution<double> mRndDistribution;
    };
}
