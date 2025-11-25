#include <vislab/field/iunsteady_vector_field3d.hpp>

namespace vislab
{
    IUnsteadyVectorField3d::IUnsteadyVectorField3d(const BoundingBox& domain)
        : IVectorField3d(domain)
    {
    }

    IUnsteadyVectorField3d::Value IUnsteadyVectorField3d::sample_dt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 0, 0, 1)); }
    IUnsteadyVectorField3d::Value IUnsteadyVectorField3d::sample_dxt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(1, 0, 0, 1)); }
    IUnsteadyVectorField3d::Value IUnsteadyVectorField3d::sample_dyt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 1, 0, 1)); }
    IUnsteadyVectorField3d::Value IUnsteadyVectorField3d::sample_dzt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 0, 1, 1)); }
    IUnsteadyVectorField3d::Value IUnsteadyVectorField3d::sample_dtt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 0, 0, 2)); }
    IUnsteadyVectorField3d::Value IUnsteadyVectorField3d::sample_dxxt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(2, 0, 0, 1)); }
    IUnsteadyVectorField3d::Value IUnsteadyVectorField3d::sample_dxyt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(1, 1, 0, 1)); }
    IUnsteadyVectorField3d::Value IUnsteadyVectorField3d::sample_dxzt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(1, 0, 1, 1)); }
    IUnsteadyVectorField3d::Value IUnsteadyVectorField3d::sample_dyyt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 2, 0, 1)); }
    IUnsteadyVectorField3d::Value IUnsteadyVectorField3d::sample_dyzt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 1, 1, 1)); }
    IUnsteadyVectorField3d::Value IUnsteadyVectorField3d::sample_dzzt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 0, 2, 1)); }
    IUnsteadyVectorField3d::Value IUnsteadyVectorField3d::sample_dxtt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(1, 0, 0, 2)); }
    IUnsteadyVectorField3d::Value IUnsteadyVectorField3d::sample_dytt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 1, 0, 2)); }
    IUnsteadyVectorField3d::Value IUnsteadyVectorField3d::sample_dztt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 0, 1, 2)); }
    IUnsteadyVectorField3d::Value IUnsteadyVectorField3d::sample_dttt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 0, 0, 3)); }

    IUnsteadyVectorField3d::Matrix IUnsteadyVectorField3d::sampleJacobian_dt(const DomainCoord& coord) const
    {
        Vector dx = sample_dxt(coord);
        Vector dy = sample_dyt(coord);
        Vector dz = sample_dzt(coord);
        Matrix result;
        result << dx, dy, dz;
        return result;
    }

    IUnsteadyVectorField3d::Value IUnsteadyVectorField3d::sampleRateOfAcceleration(const DomainCoord& coord) const
    {
        Matrix J   = sampleJacobian(coord);
        Matrix Jx  = sampleJacobian_dx(coord);
        Matrix Jy  = sampleJacobian_dy(coord);
        Matrix Jz  = sampleJacobian_dz(coord);
        Matrix Jt  = sampleJacobian_dt(coord);
        Vector v   = sample(coord);
        Vector vt  = sample_dt(coord);
        Vector vtt = sample_dtt(coord);
        return (Jx * v.x() + Jy * v.y() + Jz * v.z() + J * J + Jt) * v + J * vt + Jt * v + vtt;
    }

    IUnsteadyVectorField3d::Value IUnsteadyVectorField3d::sampleFeatureFlowField(const DomainCoord& coord) const
    {
        Matrix J  = sampleJacobian(coord);
        Vector vt = sample_dt(coord);
        return -J.inverse() * vt;
    }
}
