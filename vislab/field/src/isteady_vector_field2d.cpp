#include <vislab/field/isteady_vector_field2d.hpp>

namespace vislab
{
    ISteadyVectorField2d::ISteadyVectorField2d(const BoundingBox& domain)
        : IVectorField2d(domain)
    {
    }

    ISteadyVectorField2d::Value ISteadyVectorField2d::sampleRateOfAcceleration(const DomainCoord& coord) const
    {
        Matrix J  = sampleJacobian(coord);
        Matrix Jx = sampleJacobian_dx(coord);
        Matrix Jy = sampleJacobian_dy(coord);
        Vector v  = sample(coord);
        return (Jx * v.x() + Jy * v.y() + J * J) * v;
    }
}
