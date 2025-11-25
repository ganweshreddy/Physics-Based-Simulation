#pragma once

#include "isurface.hpp"

namespace vislab
{
    /**
     * @brief Basic interface for point geometry.
     */
    class ISurfaces : public Data
    {
        VISLAB_INTERFACE(ISurfaces, Data)

    public:
        /**
         * @brief Gets the number of surfaces in this container.
         * @return Number of surfaces.
         */
        virtual size_t getNumberOfSurfaces() const = 0;

        /**
         * @brief Deletes all surfaces (calls destructors).
         */
        virtual void clear() = 0;

        /**
         * @brief Recomputes the bounding box from the vertex buffers.
         */
        virtual void recomputeBoundingBox() = 0;

        /**
         * @brief Gets a surface.
         * @param index Index of surface to get.
         * @return Surface to get.
         */
        std::shared_ptr<ISurface> getSurface(size_t index);

        /**
         * @brief Gets a surface.
         * @param index Index of surface to get.
         * @return Surface to get.
         */
        std::shared_ptr<const ISurface> getSurface(size_t index) const;

        /**
         * @brief Adds a surface that matches the type of the derived class.
         * @param surface Surface to add.
         */
        void addSurface(std::shared_ptr<ISurface> surface);

        /**
         * @brief Tests if two surface geometries are equal.
         * @param other Surface to compare with.
         * @return True, if surfaces are equal.
         */
        virtual bool isEqual(const ISurfaces* other) const = 0;

    protected:
        /**
         * @brief Gets a surface.
         * @param index Index of surface to get.
         * @return Surface to get.
         */
        virtual std::shared_ptr<ISurface> getSurfaceImpl(size_t index) = 0;

        /**
         * @brief Gets a surface.
         * @param index Index of surface to get.
         * @return Surface to get.
         */
        virtual std::shared_ptr<const ISurface> getSurfaceImpl(size_t index) const = 0;

        /**
         * @brief Adds a surface that matches the type of the derived class.
         * @param surface Surface to add.
         */
        virtual void addSurfaceImpl(std::shared_ptr<ISurface> surface) = 0;
    };
}
