#pragma once

#include "surfaces_fwd.hpp"

#include "base_surfaces.hpp"
#include "surface.hpp"

#include <vislab/core/array.hpp>
#include <vislab/core/traits.hpp>

namespace vislab
{
    /**
     * @brief Set of surface geometry with specific internal type.
     * @tparam TArrayType Internal array type used to store vertices.
     */
    template <typename TArrayType>
    class Surfaces : public BaseSurfaces<TArrayType::Dimensions>
    {
        VISLAB_OBJECT(Surfaces, BaseSurfaces<TArrayType::Dimensions>)

    public:
        /**
         * @brief Type of the bounding box.
         */
        using BoundingBox = Eigen::AlignedBox<double, TArrayType::Dimensions>;

        /**
         * @brief Number of dimensions for vertex positions.
         */
        static const int Dimensions = TArrayType::Dimensions;

        /**
         * @brief Internal array type used to store the vertices.
         */
        using ArrayType = TArrayType;

        /**
         * @brief Internal surface type.
         */
        using SurfaceType = Surface<TArrayType>;

        /**
         * @brief Constructor.
         */
        Surfaces()
        {
        }

        /**
         * @brief Copy-constructor.
         * @param other Geometry to copy from.
         */
        Surfaces(const Surfaces& other)
            : BaseSurfaces<TArrayType::Dimensions>(other)
        {
            for (auto surface : other.surfaces)
                surfaces.push_back(std::shared_ptr<SurfaceType>(surface->clone()));
        }

        /**
         * @brief Destructor.
         */
        virtual ~Surfaces() { clear(); }

        /**
         * @brief Creates a new surface and returns a pointer to it. This function is thread-safe.
         * @return Created surface.
         */
        std::shared_ptr<SurfaceType> createSurface()
        {
            std::shared_ptr<SurfaceType> surface;
#ifdef NDEBUG
#pragma omp critical
#endif
            {
                surfaces.emplace_back(std::make_shared<SurfaceType>());
                surface = surfaces.back();
            }
            return surface;
        }

        /**
         * @brief Creates a new surface, which is copied from an existing surface and returns a pointer to it. This function is thread-safe.
         * @param copy Surface to copy.
         * @return Created surface.
         */
        std::shared_ptr<SurfaceType> createSurface(const SurfaceType& copy)
        {
            std::shared_ptr<SurfaceType> surface;
#ifdef NDEBUG
#pragma omp critical
#endif
            {
                surfaces.emplace_back(std::make_shared<SurfaceType>(copy));
                surface = surfaces.back();
            }
            return surface;
        }

        /**
         * @brief Gets the number of surfaces in this container.
         * @return Number of surfaces.
         */
        size_t getNumberOfSurfaces() const override { return surfaces.size(); }

        /**
         * @brief Deletes all surfaces
         */
        void clear() override
        {
            surfaces.clear();
        }

        /**
         * @brief Gets a surface.
         * @param index Index of the surface to get.
         * @return Surface to get.
         */
        std::shared_ptr<SurfaceType> getSurface(size_t index)
        {
            return std::dynamic_pointer_cast<SurfaceType>(getSurfaceImpl(index));
        }

        /**
         * @brief Gets a surface.
         * @param index Index of the surface to get.
         * @return Surface to get.
         */
        std::shared_ptr<const SurfaceType> getSurface(size_t index) const
        {
            return std::dynamic_pointer_cast<const SurfaceType>(getSurfaceImpl(index));
        }

        /**
         * @brief Adds a surface that matches the type of the derived class.
         * @param surface Surface to add
         */
        void addSurface(std::shared_ptr<SurfaceType> surface)
        {
#ifdef NDEBUG
#pragma omp critical
#endif
            {
                surfaces.push_back(surface);
            }
        }

        /**
         * @brief Collection of all the surfaces in the surface set
         */
        std::vector<std::shared_ptr<SurfaceType>> surfaces;

        /**
         * @brief Tests if two surface geometries are equal.
         * @param other Surface geometry to compare with.
         * @return True, if the surfaces are equal.
         */
        bool isEqual(const ISurfaces* other) const override
        {
            if (!BaseSurfaces<TArrayType::Dimensions>::isEqual(other))
                return false;
            auto otherTyped = dynamic_cast<const Surfaces*>(other);
            if (otherTyped == nullptr)
                return false;
            if (surfaces.size() != otherTyped->surfaces.size())
                return false;
            for (size_t i = 0; i < surfaces.size(); ++i)
                if (!surfaces[i]->isEqual(otherTyped->surfaces[i].get()))
                    return false;
            return true;
        }

    private:
        /**
         * @brief Gets a surface.
         * @param index Index of the surface to get.
         * @return Surface to get.
         */
        std::shared_ptr<ISurface> getSurfaceImpl(size_t index) override { return surfaces[index]; }

        /**
         * @brief Gets a surface.
         * @param index Index of the surface to get.
         * @return Surface to get.
         */
        std::shared_ptr<const ISurface> getSurfaceImpl(size_t index) const override { return surfaces[index]; }

        /**
         * @brief Adds a surface that matches the type of the derived class.
         * @param surface Surface to add.
         */
        void addSurfaceImpl(std::shared_ptr<ISurface> surface) override
        {
            auto l = std::dynamic_pointer_cast<SurfaceType>(surface);
            if (l != nullptr)
                surfaces.push_back(l);
        }
    };
}
