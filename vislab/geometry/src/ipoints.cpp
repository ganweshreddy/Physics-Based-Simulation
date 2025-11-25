#include <vislab/geometry/attributes.hpp>
#include <vislab/geometry/ipoints.hpp>

namespace vislab
{
    IPoints::IPoints()
        : mAttributes(std::make_shared<Attributes>())
    {
    }

    IPoints::IPoints(const IPoints& other)
        : mAttributes(other.mAttributes->clone())
    {
    }
    IPoints::~IPoints() {}

    std::shared_ptr<IArray> IPoints::getVertices() { return getVerticesImpl(); }
    std::shared_ptr<const IArray> IPoints::getVertices() const { return getVerticesImpl(); }

    Eigen::Index IPoints::getNumVertices() const { return getVertices()->getSize(); }
    const Eigen::Index IPoints::getDimensions() const { return getVertices()->getNumberOfComponents(); }
    size_t IPoints::getNumAttributes() const { return mAttributes->getSize(); }

    std::shared_ptr<Attributes> IPoints::getAttributes() { return mAttributes; }
    std::shared_ptr<const Attributes> IPoints::getAttributes() const { return mAttributes; }
    void IPoints::setAttributes(std::shared_ptr<Attributes> attributes) { mAttributes = attributes; }

    void IPoints::clear()
    {
        auto vertices = getVertices();
        vertices->clear();
        mAttributes->clear();
    }

    void IPoints::append(const IPoints* points)
    {
        auto vertices = getVertices();
        vertices->append(points->getVertices().get());
        assert(this->getNumAttributes() == points->getNumAttributes());
        for (size_t a = 0; a < this->getNumAttributes(); ++a)
            mAttributes->getByIndex(a)->append(points->getAttributes()->getByIndex(a).get());
    }

    void IPoints::prepend(const IPoints* points)
    {
        auto vertices = getVertices();
        vertices->prepend(points->getVertices().get());
        assert(this->getNumAttributes() == points->getNumAttributes());
        for (size_t a = 0; a < this->getNumAttributes(); ++a)
            mAttributes->getByIndex(a)->prepend(points->getAttributes()->getByIndex(a).get());
    }

    void IPoints::removeFirst(size_t n)
    {
        auto vertices = getVertices();
        vertices->removeFirst(n);
        for (size_t a = 0; a < this->getNumAttributes(); ++a)
            mAttributes->getByIndex(a)->removeFirst(n);
    }

    void IPoints::removeLast(size_t n)
    {
        auto vertices = getVertices();
        vertices->removeLast(n);
        for (size_t a = 0; a < this->getNumAttributes(); ++a)
            mAttributes->getByIndex(a)->removeLast(n);
    }

    void IPoints::reverse()
    {
        auto vertices = getVertices();
        vertices->reverse();
        for (size_t a = 0; a < this->getNumAttributes(); ++a)
            mAttributes->getByIndex(a)->reverse();
    }

    bool IPoints::isEqual(const IPoints* other) const
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
}
