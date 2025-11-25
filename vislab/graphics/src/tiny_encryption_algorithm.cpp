#include <vislab/graphics/tiny_encryption_algorithm.hpp>

namespace vislab
{
    uint32_t TinyEncryptionAlgorithm::sample_tea_32(uint32_t v0, uint32_t v1, int rounds)
    {
        uint32_t sum = 0;
        for (int i = 0; i < rounds; ++i)
        {
            sum += 0x9e3779b9;
            v0 += ((v1 << 4) + 0xa341316c) ^ (v1 + sum) ^ ((v1 >> 5) + 0xc8013ea4);
            v1 += ((v0 << 4) + 0xad90777d) ^ (v0 + sum) ^ ((v0 >> 5) + 0x7e95761e);
        }
        return v1;
    }

    uint64_t TinyEncryptionAlgorithm::sample_tea_64(uint32_t v0, uint32_t v1, int rounds)
    {
        uint32_t sum = 0;
        for (int i = 0; i < rounds; ++i)
        {
            sum += 0x9e3779b9;
            v0 += ((v1 << 4) + 0xa341316c) ^ (v1 + sum) ^ ((v1 >> 5) + 0xc8013ea4);
            v1 += ((v0 << 4) + 0xad90777d) ^ (v0 + sum) ^ ((v0 >> 5) + 0x7e95761e);
        }
        return uint64_t(v0) + (uint64_t(v1) << 32);
    }

    float TinyEncryptionAlgorithm::sample_tea_float32(uint32_t v0, uint32_t v1, int rounds)
    {
        uint32_t r = (sample_tea_32(v0, v1, rounds) >> 9) | 0x3f800000u;
        return reinterpret_cast<float&>(r) - 1.f;
    }

    double TinyEncryptionAlgorithm::sample_tea_float64(uint32_t v0, uint32_t v1, int rounds)
    {
        uint64_t r = (sample_tea_64(v0, v1, rounds) >> 12) | 0x3ff0000000000000ull;
        return reinterpret_cast<double&>(r) - 1.0;
    }
}
