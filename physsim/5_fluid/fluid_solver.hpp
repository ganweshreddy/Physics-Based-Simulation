#include <vislab/core/array_fwd.hpp>
#include <vislab/field/regular_field_fwd.hpp>

#include <Eigen/Eigen>
#include <memory>
#include <set>

namespace vislab
{
    class ColormapTexture;
}

namespace physsim
{
    /**
     * @brief Grid-based fluid solver for 2D domains.
     */
    struct FluidSolver
    {
    public:
        /**
         * @brief Enumeration of pressure solvers
         */
        enum class EPressureSolver
        {
            Iterative,
            Linear
        };

        /**
         * @brief Enumeration of boundary settings.
         */
        enum class EBoundary
        {
            Open,
            Closed
        };

        /**
         * @brief Constructor.
         * @param domain Bounding box of the domain.
         * @param resolution Number of cells for each dimension.
         */
        FluidSolver(const Eigen::AlignedBox2d& domain, const Eigen::Vector2i& resolution);

        /**
         * @brief Advance to next time step.
         */
        void advance();

        /**
         * @brief Resets the fluid flow back into the initial state.
         */
        void reset();

        /**
         * @brief Selection of the pressure solver to use.
         */
        EPressureSolver pressureSolver;

        /**
         * @brief Setting for the boundary.
         */
        EBoundary boundary;

        /**
         * @brief Accuracy for the iterative pressure solver.
         */
        float accuracy;

        /**
         * @brief Maximum number of iterations for the iterative pressure solver.
         */
        int iterations;

        /**
         * @brief Numerical integration step size.
         */
        float stepSize;

        /**
         * @brief Sets the region in which smoke is introduced each frame.
         */
        Eigen::AlignedBox2d sourceRegion;

        /**
         * @brief Texture that contains the density scalar field after color mapping.
         */
        std::shared_ptr<vislab::ColormapTexture> densityTexture;

    private:
        /**
         * @brief Sets smoke density in the source region.
         */
        void
        applyDensitySource();

        /**
         * @brief Add a vertical buoyancy force.
         */
        void addBuoyancy();

        /**
         * @brief Add the horizontal and vertical acceleration to the velocity field
         */
        void addAcceleration();

        /**
         * @brief Advects the flow into the next time step and applies boundary conditions on the vector field.
         */
        void solveAdvection();

        /**
         * @brief Perform pressure correction to make flow incompressible.
         */
        void solvePressure();

        /**
         * @brief Open domain: sets the normal component of the velocity on the outer wall, such that the numerical derivative of the velocity component at the pressure location is zero, i.e., the velocity field continues outwards smoothly.
         */
        void setNormalNeumann();

        /**
         * @brief Bounded domain: boundary condition that sets the normal component of the velocity on the outer wall to zero.
         */
        void setNormalDirichlet();

        /**
         * @brief Boundary condition that sets tangential component of the velocity to zero. This is the no-slip condition.
         */
        void setTangentialNoSlip();

        /**
         * @brief Computes the divergence on each grid point.
         */
        void computeDivergence();

        /**
         * @brief Implementation of iterative pressure solver.
         */
        void solvePoissonIterative();

        /**
         * @brief Implementation of linear pressure solver.
         */
        void solvePoissonLinear();

        /**
         * @brief Performs the pressure projection.
         */
        void correctVelocity();

        /**
         * @brief Advects velocities to the next time step.
         */
        void advectVelocity();

        /**
         * @brief Advects densities to the next time step.
         */
        void advectDensity();

        /**
         * @brief Builds Laplace operator
         * @param A Sparse output matrix.
         */
        void buildLaplace2d(Eigen::SparseMatrix<double>& A);

        /**
         * @brief Domain bounding box.
         */
        Eigen::AlignedBox2d mDomain;

        /**
         * @brief Voxel size.
         */
        Eigen::Vector2d mSpacing;

        /**
         * @brief Number of cells per dimensions.
         */
        Eigen::Vector2i mResolution;

        /**
         * @brief Densities at voxel center.
         */
        std::shared_ptr<vislab::RegularSteadyScalarField2f> mDensity;

        /**
         * @brief Pressure at voxel center.
         */
        std::shared_ptr<vislab::RegularSteadyScalarField2f> mPressure;

        /**
         * @brief Divergence at voxel center.
         */
        std::shared_ptr<vislab::RegularSteadyScalarField2f> mDivergence;

        /**
         * @brief Horizontal velocity on the vertical grid walls.
         */
        std::shared_ptr<vislab::RegularSteadyScalarField2f> mVelocity_u;

        /**
         * @brief Vertical velocity on the horizontal grid walls.
         */
        std::shared_ptr<vislab::RegularSteadyScalarField2f> mVelocity_v;

        /**
         * @brief Vertical acceleration on the horizontal grid walls.
         */
        std::shared_ptr<vislab::RegularSteadyScalarField2f> mAcceleration_u;

        /**
         * @brief Vertical acceleration on the horizontal grid walls.
         */
        std::shared_ptr<vislab::RegularSteadyScalarField2f> mAcceleration_v;

        /**
         * @brief Densities at voxel center.
         */
        std::shared_ptr<vislab::RegularSteadyScalarField2f> mDensity_tmp;

        /**
         * @brief Horizontal velocity on the vertical grid walls.
         */
        std::shared_ptr<vislab::RegularSteadyScalarField2f> mVelocity_u_tmp;

        /**
         * @brief Vertical velocity on the horizontal grid walls.
         */
        std::shared_ptr<vislab::RegularSteadyScalarField2f> mVelocity_v_tmp;

        /**
         * @brief Per-factorization of the pressure system.
         */
        Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> mPoissonFactorization;
    };
}
