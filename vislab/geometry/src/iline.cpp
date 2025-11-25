#include <vislab/geometry/attributes.hpp>
#include <vislab/geometry/iline.hpp>

namespace vislab
{
    ILine::ILine()
        : mAttributes(std::make_shared<Attributes>())
    {
    }

    ILine::ILine(const ILine& other)
        : mAttributes(std::shared_ptr<Attributes>(other.mAttributes->clone()))
    {
    }
    ILine::~ILine() {}

    std::shared_ptr<IArray> ILine::getVertices() { return getVerticesImpl(); }
    std::shared_ptr<const IArray> ILine::getVertices() const { return getVerticesImpl(); }

    Eigen::Index ILine::getNumVertices() const { return getVertices()->getSize(); }
    const Eigen::Index ILine::getDimensions() const { return getVertices()->getNumberOfComponents(); }
    size_t ILine::getNumAttributes() const { return mAttributes->getSize(); }

    std::shared_ptr<Attributes> ILine::getAttributes() { return mAttributes; }
    std::shared_ptr<const Attributes> ILine::getAttributes() const { return mAttributes; }
    void ILine::setAttributes(std::shared_ptr<Attributes> attributes) { mAttributes = attributes; }

    void ILine::clear()
    {
        auto vertices = getVertices();
        vertices->clear();
        mAttributes->clear();
    }

    void ILine::append(const ILine* line)
    {
        auto vertices = getVertices();
        vertices->append(line->getVertices().get());
        assert(this->getNumAttributes() == line->getNumAttributes());
        for (size_t a = 0; a < this->getNumAttributes(); ++a)
            mAttributes->getByIndex(a)->append(line->getAttributes()->getByIndex(a).get());
    }

    void ILine::prepend(const ILine* line)
    {
        auto vertices = getVertices();
        vertices->prepend(line->getVertices().get());
        assert(this->getNumAttributes() == line->getNumAttributes());
        for (size_t a = 0; a < this->getNumAttributes(); ++a)
            mAttributes->getByIndex(a)->prepend(line->getAttributes()->getByIndex(a).get());
    }

    void ILine::removeFirst(size_t n)
    {
        auto vertices = getVertices();
        vertices->removeFirst(n);
        for (size_t a = 0; a < this->getNumAttributes(); ++a)
            mAttributes->getByIndex(a)->removeFirst(n);
    }

    void ILine::removeLast(size_t n)
    {
        auto vertices = getVertices();
        vertices->removeLast(n);
        for (size_t a = 0; a < this->getNumAttributes(); ++a)
            mAttributes->getByIndex(a)->removeLast(n);
    }

    void ILine::reverse()
    {
        auto vertices = getVertices();
        vertices->reverse();
        for (size_t a = 0; a < this->getNumAttributes(); ++a)
            mAttributes->getByIndex(a)->reverse();
    }

    bool ILine::isEqual(const ILine* other) const
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
