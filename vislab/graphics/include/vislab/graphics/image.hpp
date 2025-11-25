#pragma once

#include "image_fwd.hpp"

#include "base_image.hpp"

#include <vislab/core/array.hpp>
#include <vislab/core/traits.hpp>

#include <memory>

namespace vislab
{
    /**
     * @brief Class that stores 2D image data.
     * @tparam TArrayType Internal array type used to store the pixel values.
     */
    template <typename TArrayType>
    class Image : public BaseImage<TArrayType::Dimensions>
    {
        VISLAB_OBJECT(Image, BaseImage<TArrayType::Dimensions>)

    public:
        /**
         * @brief Internal array type used to store the pixels.
         */
        using ArrayType = TArrayType;

        /**
         * @brief Data type for a single pixel.
         */
        using Element = typename TArrayType::Element;

        /**
         * @brief Empty default constructor.
         */
        Image()
            : mResolution(Eigen::Vector2i(0, 0))
            , mData(nullptr)
        {
        }

        /**
         * @brief Constructor of image with a certain size.
         * @param resolution Image resolution.
         */
        Image(const Eigen::Vector2i& resolution)
            : mResolution(resolution)
        {
            mData = std::make_shared<TArrayType>();
            mData->setSize(resolution.x() * resolution.y());
        }

        /**
         * @brief Constructor of image with a certain size.
         * @param width Width of image.
         * @param height Height of image.
         */
        Image(const int& width, const int& height)
            : mResolution(Eigen::Vector2i(width, height))
        {
            mData = std::make_shared<TArrayType>();
            mData->setSize(width * height);
        }

        /**
         * @brief Destructor.
         */
        ~Image() {}

        /**
         * @brief Sets the value at a certain pixel.
         * @param px x-coordinate of pixel.
         * @param py y-coordinate of pixel.
         * @param value Value to set.
         */
        void setValue(const int& px, const int& py, const Element& value) { setValue((Eigen::Index)py * mResolution.x() + px, value); }

        /**
         * @brief Sets the value at a certain pixel.
         * @param pixel Pixel index.
         * @param value Value to set.
         */
        void setValue(const Eigen::Vector2i& pixel, const Element& value) { setValue(pixel.y() * mResolution.x() + pixel.x(), value); }

        /**
         * @brief Sets the value at a certain pixel.
         * @param linearIndex Linear array index.
         * @param value Value to set.
         */
        void setValue(const Eigen::Index& linearIndex, const Element& value) { mData->setValue(linearIndex, value); }

        /**
         * @brief Gets the value at a certain pixel.
         * @param px x-coordinate of pixel.
         * @param py y-coordinate of pixel.
         * @return Value at pixel.
         */
        typename TArrayType::ConstElementRef getValue(const int& px, const int& py) const { return getValue((Eigen::Index)py * mResolution.x() + px); }

        /**
         * @brief Gets the value at a certain pixel.
         * @param pixel Pixel index.
         * @return Value at pixel.
         */
        typename TArrayType::ConstElementRef getValue(const Eigen::Vector2i& pixel) const { return getValue(pixel.y() * mResolution.x() + pixel.x()); }

        /**
         * @brief Gets the value at a certain pixel.
         * @param linearIndex Linear array index.
         * @return Value at pixel.
         */
        typename TArrayType::ConstElementRef getValue(const Eigen::Index& linearIndex) const { return mData->getValue(linearIndex); }

        /**
         * @brief Gets the physical coordinate of a grid point in double precision, identified by its linear index.
         * @param linearIndex Linear array index.
         * @return Value at pixel.
         */
        Eigen::Vector<double, TArrayType::Dimensions> getValueDouble(const Eigen::Index& linearIndex) const override { return mData->getValueDouble(linearIndex); }

        /**
         * @brief Sets the value at a certain pixel in double precision.
         * @param linearIndex Linear array index.
         * @param value Value to set.
         */
        void setValueDouble(const Eigen::Index& linearIndex, const Eigen::Vector<double, TArrayType::Dimensions>& value) override { mData->setValueDouble(linearIndex, value); }

        /**
         * @brief Gets the resolution of the image.
         * @return Image resolution.
         */
        const Eigen::Vector2i& getResolution() const override { return mResolution; }

        /**
         * @brief Sets the resolution of the image. The internal linear image data array gets resized.
         * @param resolution New image resolution.
         */
        void setResolution(const Eigen::Vector2i& resolution) override
        {
            mResolution = resolution;
            if (mData == nullptr)
                mData = std::make_shared<TArrayType>();
            mData->setSize(resolution.x() * resolution.y());
        }

        /**
         * @brief Sets the resolution of the image.
         * @param resx New width of the image.
         * @param resy New height of the image.
         */
        void setResolution(const int resx, const int resy) override
        {
            setResolution(Eigen::Vector2i(resx, resy));
        }

        /**
         * @brief Gets the total number of pixels.
         * @return Number of pixels.
         */
        Eigen::Index getNumPixels() const override { return mResolution.x() * mResolution.y(); }

        /**
         * @brief Sets all values to zero.
         */
        void setZero() override
        {
            if (mData != nullptr)
                mData->setZero();
        }

        /**
         * @brief Gets the array that stores the data in the image.
         * @return Internal data array.
         */
        std::shared_ptr<TArrayType> getArray() { return mData; }

        /**
         * @brief Gets the array that stores the data in the image.
         * @return Internal data array.
         */
        std::shared_ptr<const TArrayType> getArray() const { return mData; }

    private:
        /**
         * @brief Image resolution.
         */
        Eigen::Vector2i mResolution;

        /**
         * @brief Data stored in the image.
         */
        std::shared_ptr<TArrayType> mData;
    };
}
