#include <Eigen/Eigen>

namespace vislab
{
    /**
     * @brief Class that represents a ray.
     */
    template <typename TVector>
    class Ray;

    /**
     * @brief Ray in two-dimensional space.
     */
    using Ray2d = Ray<Eigen::Vector2d>;

    /**
     * @brief Ray in three-dimensional space.
     */
    using Ray3d = Ray<Eigen::Vector3d>;
}
