#include <vislab/field/iunsteady_vector_field2d.hpp>

namespace vislab
{
    IUnsteadyVectorField2d::IUnsteadyVectorField2d(const BoundingBox& domain)
        : IVectorField2d(domain)
    {
    }

    IUnsteadyVectorField2d::Value IUnsteadyVectorField2d::sample_dt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 0, 1)); }
    IUnsteadyVectorField2d::Value IUnsteadyVectorField2d::sample_dxt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(1, 0, 1)); }
    IUnsteadyVectorField2d::Value IUnsteadyVectorField2d::sample_dyt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 1, 1)); }
    IUnsteadyVectorField2d::Value IUnsteadyVectorField2d::sample_dtt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 0, 2)); }
    IUnsteadyVectorField2d::Value IUnsteadyVectorField2d::sample_dxxt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(2, 0, 1)); }
    IUnsteadyVectorField2d::Value IUnsteadyVectorField2d::sample_dxyt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(1, 1, 1)); }
    IUnsteadyVectorField2d::Value IUnsteadyVectorField2d::sample_dyyt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 2, 1)); }
    IUnsteadyVectorField2d::Value IUnsteadyVectorField2d::sample_dxtt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(1, 0, 2)); }
    IUnsteadyVectorField2d::Value IUnsteadyVectorField2d::sample_dytt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 1, 2)); }
    IUnsteadyVectorField2d::Value IUnsteadyVectorField2d::sample_dttt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 0, 3)); }

    IUnsteadyVectorField2d::Matrix IUnsteadyVectorField2d::sampleJacobian_dt(const DomainCoord& coord) const
    {
        Vector dx = sample_dxt(coord);
        Vector dy = sample_dyt(coord);
        Matrix result;
        result << dx, dy;
        return result;
    }

    IUnsteadyVectorField2d::Matrix IUnsteadyVectorField2d::sampleMaterialDerivativeJacobian(const DomainCoord& coord) const
    {
        Vector v  = sample(coord);
        Matrix Jx = sampleJacobian_dx(coord);
        Matrix Jy = sampleJacobian_dy(coord);
        Matrix Jt = sampleJacobian_dt(coord);
        return Jx * v.x() + Jy * v.y() + Jt;
    }

    IUnsteadyVectorField2d::Value IUnsteadyVectorField2d::sampleRateOfAcceleration(const DomainCoord& coord) const
    {
        Matrix J   = sampleJacobian(coord);
        Matrix Jx  = sampleJacobian_dx(coord);
        Matrix Jy  = sampleJacobian_dy(coord);
        Matrix Jt  = sampleJacobian_dt(coord);
        Vector v   = sample(coord);
        Vector vt  = sample_dt(coord);
        Vector vtt = sample_dtt(coord);
        return (Jx * v.x() + Jy * v.y() + J * J + Jt) * v + J * vt + Jt * v + vtt;
    }

    IUnsteadyVectorField2d::Value IUnsteadyVectorField2d::sampleFeatureFlowField(const DomainCoord& coord) const
    {
        Matrix J  = sampleJacobian(coord);
        Vector vt = sample_dt(coord);
        return -J.completeOrthogonalDecomposition().solve(vt);
    }
}
