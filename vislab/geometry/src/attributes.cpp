#include <vislab/geometry/attributes.hpp>

#include <vislab/core/iarray.hpp>

namespace vislab
{
    Attributes::Attributes() {}

    Attributes::Attributes(const Attributes& other)
    {
        mAttributes.resize(other.getSize());
        for (size_t i = 0; i < mAttributes.size(); ++i)
        {
            mAttributes[i] = std::shared_ptr<IArray>(other.getByIndex(i)->clone());
        }
    }
    Attributes::~Attributes() {}

    size_t Attributes::getSize() const { return mAttributes.size(); }

    std::shared_ptr<IArray> Attributes::getByIndex(const size_t& index) { return mAttributes[index]; }
    std::shared_ptr<const IArray> Attributes::getByIndex(const size_t& index) const { return mAttributes[index]; }

    std::shared_ptr<IArray> Attributes::getByName(const std::string& name)
    {
        for (auto& attr : mAttributes)
            if (attr->getName() == name)
                return attr;
        return nullptr;
    }

    std::shared_ptr<const IArray> Attributes::getByName(const std::string& name) const
    {
        for (auto& attr : mAttributes)
            if (attr->getName() == name)
                return attr;
        return nullptr;
    }

    void Attributes::clear()
    {
        mAttributes.clear();
    }

    void Attributes::append(const Attributes* attributes)
    {
        assert(this->getSize() == attributes->getSize());
        for (size_t a = 0; a < this->getSize(); ++a)
            this->getByIndex(a)->append(attributes->getByIndex(a).get());
    }

    void Attributes::prepend(const Attributes* attributes)
    {
        assert(this->getSize() == attributes->getSize());
        for (size_t a = 0; a < this->getSize(); ++a)
            this->getByIndex(a)->prepend(attributes->getByIndex(a).get());
    }

    void Attributes::removeFirst(size_t n)
    {
        for (size_t a = 0; a < this->getSize(); ++a)
            this->getByIndex(a)->removeFirst(n);
    }

    void Attributes::removeLast(size_t n)
    {
        for (size_t a = 0; a < this->getSize(); ++a)
            this->getByIndex(a)->removeLast(n);
    }

    void Attributes::reverse()
    {
        for (size_t a = 0; a < this->getSize(); ++a)
            this->getByIndex(a)->reverse();
    }

    bool Attributes::isEqual(const Attributes* other) const
    {
        if (this->getSize() != other->getSize())
            return false;
        for (size_t i = 0; i < other->getSize(); ++i)
            if (!this->getByIndex(i)->isEqual(other->getByIndex(i).get()))
                return false;
        return true;
    }

    const std::string& Attributes::getName(const size_t& attributeIndex) const
    {
        assert(attributeIndex >= 0 && attributeIndex < mAttributes.size() && "Attribute with index does not exist.");
        return mAttributes[attributeIndex]->getName();
    }
}
