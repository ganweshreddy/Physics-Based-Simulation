#include <vislab/graphics/sampler.hpp>

namespace vislab
{
    Sampler::Sampler()
        : mSampleCount(4)
        , mSampleIndex(0)
        , mBaseSeed(0)
    {
    }

    Sampler::Sampler(const Sampler& other)
        : mSampleCount(other.mSampleCount)
        , mSampleIndex(other.mSampleIndex)
        , mBaseSeed(other.mBaseSeed)
    {
    }

    void Sampler::seed(uint64_t seed_offset)
    {
        mSampleIndex = 0;
    }

    void Sampler::advance()
    {
        assert(mSampleIndex < mSampleCount);
        mSampleIndex++;
    }

    uint32_t Sampler::getSampleIndex() const { return mSampleIndex; }
    uint32_t Sampler::getSampleCount() const { return mSampleCount; }
    uint64_t Sampler::getBaseSeed() const { return mBaseSeed; }
    void Sampler::setBaseSeed(uint64_t baseSeed) { mBaseSeed = baseSeed; }
}
