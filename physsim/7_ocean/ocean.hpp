#include <wave.hpp>

#include <vislab/core/array_fwd.hpp>

#include <memory>
#include <vector>

namespace vislab
{
    class Mesh;
}

namespace physsim
{
    /**
     * @brief Helper class for the implementation of ocean simulation based on Gerstner waves.
     */
    struct Ocean
    {
    public:
        /**
         * @brief Constructor for a new ocean.
         * @param resolution Resolution of the grid.
         * @param origin Origin of the mesh grid.
         * @param axis1 Horizontal axis for the grid (x coordinate spacing).
         * @param axis2 Vertical axis for the grid (y coordinate spacing).
         */
        Ocean(const Eigen::Vector2i& resolution, const Eigen::Vector3d& origin, const Eigen::Vector3d& axis1, const Eigen::Vector3d& axis2);

        /**
         * @brief Advance to next time step.
         * @param totalTime Total time in milliseconds since the beginning of the first frame.
         */
        void advance(double totalTime);

        /**
         * @brief Resolution of the ocean grid.
         * @return Grid resolution.
         */
        const Eigen::Vector2i& resolution() const;

        /**
         * @brief Mesh that stores the geometry.
         */
        std::shared_ptr<vislab::Mesh> mesh;

        /**
         * @brief Collection of waves to add up.
         */
        std::vector<Wave> waves;

        /**
         * @brief Evaluates the displaced position at a specific domain position.
         * @param pos Domain position of the vertex.
         * @param pos Time to evaluate in milliseconds.
         * @return Displaced location of the vertex.
         */
        Eigen::Vector3d position(const Eigen::Vector3d& pos, double t) const;

    private:
        /**
         * @brief Resolution of the grid.
         */
        Eigen::Vector2i mResolution;

        /**
         * @brief Origin of the mesh grid.
         */
        Eigen::Vector3d mOrigin;

        /**
         * @brief Horizontal axis for the grid (x coordinate spacing).
         */
        Eigen::Vector3d mAxis1;

        /**
         * @brief Vertical axis for the grid (y coordinate spacing).
         */
        Eigen::Vector3d mAxis2;
    };
}
