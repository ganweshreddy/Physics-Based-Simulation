#include <vislab/field/iunsteady_scalar_field3d.hpp>

namespace vislab
{
    IUnsteadyScalarField3d::IUnsteadyScalarField3d(const BoundingBox& domain)
        : IScalarField3d(domain)
    {
    }

    IUnsteadyScalarField3d::Value IUnsteadyScalarField3d::sample_dt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 0, 0, 1)); }
    IUnsteadyScalarField3d::Value IUnsteadyScalarField3d::sample_dxt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(1, 0, 0, 1)); }
    IUnsteadyScalarField3d::Value IUnsteadyScalarField3d::sample_dyt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 1, 0, 1)); }
    IUnsteadyScalarField3d::Value IUnsteadyScalarField3d::sample_dzt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 0, 1, 1)); }
    IUnsteadyScalarField3d::Value IUnsteadyScalarField3d::sample_dtt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 0, 0, 2)); }
    IUnsteadyScalarField3d::Value IUnsteadyScalarField3d::sample_dxxt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(2, 0, 0, 1)); }
    IUnsteadyScalarField3d::Value IUnsteadyScalarField3d::sample_dxyt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(1, 1, 0, 1)); }
    IUnsteadyScalarField3d::Value IUnsteadyScalarField3d::sample_dxzt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(1, 0, 1, 1)); }
    IUnsteadyScalarField3d::Value IUnsteadyScalarField3d::sample_dyyt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 2, 0, 1)); }
    IUnsteadyScalarField3d::Value IUnsteadyScalarField3d::sample_dyzt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 1, 1, 1)); }
    IUnsteadyScalarField3d::Value IUnsteadyScalarField3d::sample_dzzt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 0, 2, 1)); }
    IUnsteadyScalarField3d::Value IUnsteadyScalarField3d::sample_dxtt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(1, 0, 0, 2)); }
    IUnsteadyScalarField3d::Value IUnsteadyScalarField3d::sample_dytt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 1, 0, 2)); }
    IUnsteadyScalarField3d::Value IUnsteadyScalarField3d::sample_dztt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 0, 1, 2)); }
    IUnsteadyScalarField3d::Value IUnsteadyScalarField3d::sample_dttt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 0, 0, 3)); }
}
