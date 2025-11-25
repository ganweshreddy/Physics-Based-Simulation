#include <vislab/field/isteady_vector_field3d.hpp>

namespace vislab
{
    ISteadyVectorField3d::ISteadyVectorField3d(const BoundingBox& domain)
        : IVectorField3d(domain)
    {
    }

    ISteadyVectorField3d::Value ISteadyVectorField3d::sampleRateOfAcceleration(const DomainCoord& coord) const
    {
        Matrix J  = sampleJacobian(coord);
        Matrix Jx = sampleJacobian_dx(coord);
        Matrix Jy = sampleJacobian_dy(coord);
        Matrix Jz = sampleJacobian_dz(coord);
        Vector v  = sample(coord);
        return (Jx * v.x() + Jy * v.y() + Jz * v.z() + J * J) * v;
    }
}
