#include <vislab/graphics/independent_sampler.hpp>

namespace vislab
{
    IndependentSampler::IndependentSampler()
        : Sampler()
    {
    }
    IndependentSampler::IndependentSampler(const IndependentSampler& other)
        : Sampler(other)
    {
    }

    void IndependentSampler::seed(uint64_t seed_offset)
    {
        Sampler::seed(seed_offset);
        mRndGenerator.seed(mBaseSeed + seed_offset);
    }

    double IndependentSampler::next_1d()
    {
        return mRndDistribution(mRndGenerator);
    }

    Eigen::Vector2d IndependentSampler::next_2d()
    {
        return Eigen::Vector2d(
            next_1d(),
            next_1d());
    }
}
