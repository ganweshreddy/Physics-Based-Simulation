#pragma once

#define VISLAB_INTERFACE(CLASS, ...)                                                                                                                                                                  \
public:                                                                                                                                                                                               \
    /**                                                                                                                                                                                               \
     * @brief Makes a deep copy of this object by using the copy-constructor of the derived class. This method is overwritten by each derived class in order to return the type of the derived class. \
     * @return Deep copy of this object.                                                                                                                                                              \
     */                                                                                                                                                                                               \
    inline CLASS* clone() const                                                                                                                                                                       \
    {                                                                                                                                                                                                 \
        return this->cloneImpl();                                                                                                                                                                     \
    }                                                                                                                                                                                                 \
                                                                                                                                                                                                      \
private:                                                                                                                                                                                              \
    /**                                                                                                                                                                                               \
     * @brief Actual implementation of the deep cloning.                                                                                                                                              \
     * @return Deep copy of this object.                                                                                                                                                              \
     */                                                                                                                                                                                               \
    virtual CLASS* cloneImpl() const = 0;

#define VISLAB_OBJECT(CLASS, ...)                                                                                                                                                                     \
public:                                                                                                                                                                                               \
    /**                                                                                                                                                                                               \
     * @brief Makes a deep copy of this object by using the copy-constructor of the derived class. This method is overwritten by each derived class in order to return the type of the derived class. \
     * @return Deep copy of this object.                                                                                                                                                              \
     */                                                                                                                                                                                               \
    inline CLASS* clone() const                                                                                                                                                                       \
    {                                                                                                                                                                                                 \
        return this->cloneImpl();                                                                                                                                                                     \
    }                                                                                                                                                                                                 \
                                                                                                                                                                                                      \
private:                                                                                                                                                                                              \
    /**                                                                                                                                                                                               \
     * @brief Actual implementation of the deep cloning.                                                                                                                                              \
     * @return Deep copy of this object.                                                                                                                                                              \
     */                                                                                                                                                                                               \
    inline virtual CLASS* cloneImpl() const                                                                                                                                                           \
    {                                                                                                                                                                                                 \
        return new CLASS(*this);                                                                                                                                                                      \
    }
