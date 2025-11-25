#include <vislab/field/iunsteady_scalar_field2d.hpp>

namespace vislab
{
    IUnsteadyScalarField2d::IUnsteadyScalarField2d(const BoundingBox& domain)
        : IScalarField2d(domain)
    {
    }

    IUnsteadyScalarField2d::Value IUnsteadyScalarField2d::sample_dt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 0, 1)); }
    IUnsteadyScalarField2d::Value IUnsteadyScalarField2d::sample_dxt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(1, 0, 1)); }
    IUnsteadyScalarField2d::Value IUnsteadyScalarField2d::sample_dyt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 1, 1)); }
    IUnsteadyScalarField2d::Value IUnsteadyScalarField2d::sample_dtt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 0, 2)); }
    IUnsteadyScalarField2d::Value IUnsteadyScalarField2d::sample_dxxt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(2, 0, 1)); }
    IUnsteadyScalarField2d::Value IUnsteadyScalarField2d::sample_dxyt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(1, 1, 1)); }
    IUnsteadyScalarField2d::Value IUnsteadyScalarField2d::sample_dyyt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 2, 1)); }
    IUnsteadyScalarField2d::Value IUnsteadyScalarField2d::sample_dxtt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(1, 0, 2)); }
    IUnsteadyScalarField2d::Value IUnsteadyScalarField2d::sample_dytt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 1, 2)); }
    IUnsteadyScalarField2d::Value IUnsteadyScalarField2d::sample_dttt(const DomainCoord& coord) const { return samplePartial(coord, DerivativeDegree(0, 0, 3)); }
}
