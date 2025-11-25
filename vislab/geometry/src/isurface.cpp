#include <vislab/geometry/isurface.hpp>

#include <vislab/geometry/attributes.hpp>

#include <vislab/core/base_array.hpp>

namespace vislab
{
    ISurface::ISurface()
        : mPrimitiveTopology(ISurface::ETopology::TriangleList)
        , mAttributes(std::make_shared<Attributes>())
    {
    }

    ISurface::ISurface(const ISurface& other)
        : mPrimitiveTopology(ISurface::ETopology::TriangleList)
        , mAttributes(std::shared_ptr<Attributes>(other.mAttributes->clone()))
    {
        if (other.mIndices)
            mIndices = std::shared_ptr<IArray1>(other.mIndices->clone());
    }

    ISurface::~ISurface() {}

    std::shared_ptr<IArray> ISurface::getVertices() { return getVerticesImpl(); }
    std::shared_ptr<const IArray> ISurface::getVertices() const { return getVerticesImpl(); }

    Eigen::Index ISurface::getNumVertices() const { return getVertices()->getSize(); }
    const Eigen::Index ISurface::getDimensions() const { return getVertices()->getNumberOfComponents(); }
    size_t ISurface::getNumAttributes() const { return mAttributes->getSize(); }

    Eigen::Index ISurface::getNumIndices() const
    {
        if (mIndices)
            return mIndices->getSize();
        else
            return 0;
    }

    std::shared_ptr<Attributes> ISurface::getAttributes() { return mAttributes; }
    std::shared_ptr<const Attributes> ISurface::getAttributes() const { return mAttributes; }
    void ISurface::setAttributes(std::shared_ptr<Attributes> attributes) { mAttributes = attributes; }

    void ISurface::clear()
    {
        auto vertices = getVertices();
        vertices->clear();
        mAttributes->clear();
    }

    bool ISurface::isEqual(const ISurface* other) const
    {
        auto vertices = getVertices();
        if (!vertices->isEqual(other->getVertices().get()))
            return false;
        if (this->getNumAttributes() != other->getNumAttributes())
            return false;
        for (size_t i = 0; i < other->getNumAttributes(); ++i)
            if (!mAttributes->getByIndex(i)->isEqual(other->getAttributes()->getByIndex(i).get()))
                return false;
        return true;
    }

    std::shared_ptr<IArray1> ISurface::getIndices() { return mIndices; }
    std::shared_ptr<const IArray1> ISurface::getIndices() const { return mIndices; }

    void ISurface::clearIndices()
    {
        mIndices = nullptr;
    }

    void ISurface::cloneIndices(const ISurface& other)
    {
        mIndices = std::shared_ptr<IArray1>(other.mIndices->clone());
    }

    const ISurface::ETopology& ISurface::getPrimitiveTopology() const { return mPrimitiveTopology; }
    void ISurface::setPrimitiveTopology(const ETopology& topology) { mPrimitiveTopology = topology; }
}
