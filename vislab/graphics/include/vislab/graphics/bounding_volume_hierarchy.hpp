#pragma once

#include "preliminary_interaction.hpp"
#include "ray.hpp"
#include "ray_box_intersect.hpp"

#include <algorithm>

namespace vislab
{
    /**
     * @brief Class that forms a bounding volume hierarchy over shapes of type TShape.
     * @tparam TShape Type of shapes stored in the hierarchy.
     * @tparam TDimensions Dimensionality of the domain.
     */
    template <class TShape, int TDimensions>
    class BoundingVolumeHierarchy
    {
    public:
        /**
         * @brief Dimensionality of the tree.
         */
        static const int Dimensions = TDimensions;

        /**
         * @brief Constructor of an empty tree.
         */
        BoundingVolumeHierarchy()
            : mRoot(NULL)
        {
            for (int i = 0; i < TDimensions; ++i)
                mSortedObjects[i] = NULL;
        }

        /**
         * @brief Copies the shapes from another tree, but the acceleration data structure is not built yet.
         */
        BoundingVolumeHierarchy(const BoundingVolumeHierarchy& other)
            : mRoot(NULL)
        {
            for (int i = 0; i < TDimensions; ++i)
                mSortedObjects[i] = NULL;
            mShapes = other.mShapes;
        }

        /**
         * @brief Destructor.
         */
        ~BoundingVolumeHierarchy()
        {
            if (mRoot != NULL)
                delete mRoot;

            int numShapes = int(mShapes.size());
            for (int i = 0; i < numShapes; i++)
            {
                for (int k = 0; k < TDimensions; ++k)
                {
                    delete mSortedObjects[k][i];
                    mSortedObjects[k][i] = NULL;
                }
            }
            for (int k = 0; k < TDimensions; ++k)
            {
                delete[] mSortedObjects[k];
                mSortedObjects[k] = NULL;
            }
            mShapes.clear();
        }

        /**
         * @brief Gets the shapes that are stored in this tree.
         * @return Vector of shapes.
         */
        const std::vector<TShape*>& getShapes() const { return mShapes; }

        /**
         * @brief Gets the shapes that are stored in this tree.
         * @return Writable vector of shapes.
         */
        std::vector<TShape*>& getShapes() { return mShapes; }

        /**
         * @brief Builds the tree from formerly added Objects
         */
        void build()
        {
            int numShapes = int(mShapes.size());
            for (int k = 0; k < TDimensions; ++k)
                mSortedObjects[k] = new Leaf*[numShapes];
            bool* isLeft = new bool[numShapes];

            // init sorted pointer lists
            for (int i = 0; i < numShapes; i++)
            {
                for (int k = 0; k < TDimensions; ++k)
                    mSortedObjects[k][i] = new Leaf(mShapes[i], &isLeft[i]);

                //mShapes[i]->recomputeBoundingBox();
            }

            for (int k = 0; k < TDimensions; ++k)
                sortObjects(mSortedObjects[k], 0, numShapes - 1, k);

            mRoot = new Node();
            Eigen::AlignedBox<double, TDimensions> rootBox;
            rootBox.setEmpty();

            for (int i = 0; i < numShapes; i++)
            {
                Leaf* leaf = mSortedObjects[0][i];
                //leaf->shape->recomputeBoundingBox();
                rootBox.extend(leaf->shape->worldBounds());
            }

            createHierarchy(mRoot, rootBox, mSortedObjects, numShapes);
            mRoot->recomputeBoundingBox();

            delete[] isLeft;
        }

        /**
         * @brief Fast ray intersection test. To obtain further information about the intersection, use ComputeSurfaceInteraction().
         * @param ray Ray to test against.
         * @return Intersection result.
         */
        PreliminaryIntersection preliminaryHit(const Ray3d& ray) const
        {
            if (mShapes.size() == 0)
                return PreliminaryIntersection();
            if (mShapes.size() == 1)
                return mShapes[0]->preliminaryHit(ray);

            if (mRoot)
                return mRoot->preliminaryHit(ray);
            else
                return PreliminaryIntersection();
        }

        /**
         * @brief Fast ray shadow test.
         * @param ray Ray to test against.
         * @return True if there was a hit.
         */
        bool anyHit(const Ray3d& ray) const
        {
            if (mShapes.size() == 0)
                return false;
            if (mShapes.size() == 1)
                return mShapes[0]->anyHit(ray);
            if (mRoot)
                return mRoot->anyHit(ray);
            else
                return false;
        }

    private:
        /**
         * @brief Shape that is inserted into the BVH tree (a leaf). Class T must have a "getBoundingBox" function.
         */
        class Leaf
        {
        public:
            /**
             * @brief Constructor.
             * @param shape Shape to be stored in this tree.
             * @param isLeft Flag that identifies whether this is a left leaf.
             */
            Leaf(TShape* shape, bool* isLeft)
                : shape(shape)
                , isLeft(isLeft)
            {
            }

            /**
             * @brief Shape that is stored in this leaf.
             */
            TShape* shape;

            /**
             * @brief Flag that determines whether this is a left leaf.
             */
            bool* isLeft;
        };

        /**
         * @brief Inner node of the bounding volume hierarchy.
         */
        class Node
        {
        public:
            /**
             * @brief Constructor that creates an empty node.
             */
            Node()
                : mLeft(NULL)
                , mRight(NULL)
            {
                mLeaf[0] = NULL;
                mLeaf[1] = NULL;
                mBounds.setEmpty();
            }

            /**
             * @brief Destructor that deletes all shapes below this node.
             */
            ~Node()
            {
                delete mLeft;
                delete mRight;
            }

            /**
             * @brief Adds a node on the left side.
             * @param child Child node to add.
             */
            void addLeft(Node* child)
            {
                if (!mLeft)
                    mLeft = child;
            }

            /**
             * @brief Adds a node on the right side.
             * @param child Child node to add.
             */
            void addRight(Node* child)
            {
                if (!mRight)
                    mRight = child;
            }

            /**
             * @brief Sets both leafs at once.
             * @param leaf1 Left leaf.
             * @param leaf2 Right leaf.
             */
            void setLeafs(const Leaf* leaf1, const Leaf* leaf2)
            {
                mLeaf[0] = leaf1;
                mLeaf[1] = leaf2;
            }

            /**
             * @brief Recomputes the bounding box from the childs and leafs.
             */
            void recomputeBoundingBox()
            {
                for (int i = 0; i < 2; ++i)
                    if (mLeaf[i])
                        mBounds.extend(mLeaf[i]->shape->worldBounds());

                if (mLeft)
                    mBounds.extend(mLeft->getBoundingBox());
                if (mRight)
                    mBounds.extend(mRight->getBoundingBox());
            }

            /**
             * @brief Gets the bounding box.
             * @return Bounding of the node.
             */
            const Eigen::AlignedBox<double, TDimensions>& getBoundingBox() const { return mBounds; }

            /**
             * @brief Fast ray intersection test. To obtain further information about the intersection, use ComputeSurfaceInteraction().
             * @param ray Ray to test against.
             * @return Intersection result.
             */
            virtual PreliminaryIntersection preliminaryHit(const Ray3d& ray) const
            {
                // Does the ray hit the bounding box?
                if (!RayBoxIntersection::anyHit(ray, mBounds))
                    return PreliminaryIntersection();

                // Make local copy of ray.
                Ray3d ray_ = ray;
                PreliminaryIntersection piLeft, piRight;

                // check left child tree
                if (mLeft != nullptr)
                    piLeft = mLeft->preliminaryHit(ray_);
                else if (mLeaf[0] != nullptr)
                    piLeft = mLeaf[0]->shape->preliminaryHit(ray_);
                if (piLeft.isValid())
                    ray_.tMax = piLeft.t;

                // check right child tree
                if (mRight != nullptr)
                    piRight = mRight->preliminaryHit(ray_);
                else if (mLeaf[1] != nullptr)
                    piRight = mLeaf[1]->shape->preliminaryHit(ray_);

                // if the right is valid, it is defintely closer than the left, since we adjusted tMax.
                if (piRight.isValid())
                    return piRight;
                return piLeft;
            }

            /**
             * @brief Checks if there is any intersection.
             * @param ray Ray to test against.
             * @return True if there was a hit.
             */
            bool anyHit(const Ray3d& ray) const
            {
                // Does the ray hit the bounding box?
                if (!RayBoxIntersection::anyHit(ray, mBounds))
                    return false;

                // hit the Object?
                if (mLeaf[0] && mLeaf[0]->shape->anyHit(ray))
                    return true;
                if (mLeaf[1] && mLeaf[1]->shape->anyHit(ray))
                    return true;

                // hit something left?
                if (mLeft && mLeft->anyHit(ray))
                    return true;

                // hit something right?
                if (mRight)
                    return mRight->anyHit(ray);

                // no hit
                return false;
            }

        private:
            /**
             * @brief Left child node.
             */
            Node* mLeft;

            /**
             * @brief Right child node.
             */
            Node* mRight;

            /**
             * @brief Bounding box of the node.
             */
            Eigen::AlignedBox<double, TDimensions> mBounds;

            /**
             * @brief The two leafs.
             */
            const Leaf* mLeaf[2];
        };

        /**
         * @brief Sorts the shapes.
         * @param sortedObjects Pointer to the shapes to sort.
         * @param left Starting index.
         * @param right Ending index.
         * @param axis Dimension over which to sort.
         */
        void sortObjects(Leaf** sortedObjects, int left, int right, int axis)
        {
            Leaf* tmpObject;
            int i, j;
            double comparisonCenter;

            if (right <= left)
                return;

            comparisonCenter = sortedObjects[left]->shape->worldBounds().center()[axis];
            i                = left + 1;
            j                = right;

            // special case of all coordinates being identical -> pick the middle
            if (sortedObjects[i]->shape->worldBounds().center()[axis] == sortedObjects[j]->shape->worldBounds().center()[axis])
            {
                i = (left + right) / 2;
                j = (left + right) / 2;
            }
            else
            {
                do
                {
                    while ((i < right) && (sortedObjects[i]->shape->worldBounds().center()[axis] <= comparisonCenter))
                        i++;

                    while (sortedObjects[j]->shape->worldBounds().center()[axis] > comparisonCenter)
                        j--;

                    if (i < j)
                    {
                        tmpObject        = sortedObjects[i];
                        sortedObjects[i] = sortedObjects[j];
                        sortedObjects[j] = tmpObject;
                    }
                } while (i < j);
            }

            if (j != left)
            {
                tmpObject           = sortedObjects[j];
                sortedObjects[j]    = sortedObjects[left];
                sortedObjects[left] = tmpObject;
            }

            sortObjects(sortedObjects, left, j - 1, axis);
            sortObjects(sortedObjects, j + 1, right, axis);
        }

        /**
         * @brief
         * @param obj
         * @param bbox
         * @param sortedObjectsList
         * @param numLeafs
         */
        void createHierarchy(Node* node, const Eigen::AlignedBox<double, TDimensions>& boundingBox, std::array<Leaf**, TDimensions> sortedObjectsList, int numLeafs)
        {
            double totalArea, minCost, splitCost, leftArea, rightArea;
            int axis, splitAxis = 0, splitIndex = 0, i, leftIndex = 0, rightIndex = 0;
            Eigen::AlignedBox<double, TDimensions>*leftBoxes, *rightBoxes;
            Eigen::AlignedBox<double, TDimensions> leftBox, rightBox;
            std::array<Leaf**, TDimensions> leftSortedObjectsList  = { NULL };
            std::array<Leaf**, TDimensions> rightSortedObjectsList = { NULL };
            Leaf** sortedObjects                                   = NULL;
            Leaf** leftSortedObjects                               = NULL;
            Leaf** rightSortedObjects                              = NULL;

            if (numLeafs <= 1)
            {
                return;
            }
            if (numLeafs == 2)
            {
                node->setLeafs(sortedObjectsList[0][0], sortedObjectsList[0][1]);
                node->recomputeBoundingBox();
                return;
            }

            leftBoxes  = new Eigen::AlignedBox<double, TDimensions>[numLeafs - 1];
            rightBoxes = new Eigen::AlignedBox<double, TDimensions>[numLeafs - 1];

            totalArea = computeArea(boundingBox);
            minCost   = std::numeric_limits<float>::max();

            for (axis = 0; axis < TDimensions; axis++)
            {
                sortedObjects = sortedObjectsList[axis];
                createLeftBoxes(leftBoxes, numLeafs, sortedObjects);
                createRightBoxes(rightBoxes, numLeafs, sortedObjects);

                for (i = 0; i <= numLeafs - 2; i++)
                {
                    leftArea  = computeArea(leftBoxes[i]);
                    rightArea = computeArea(rightBoxes[i]);
                    splitCost = (leftArea * (i + 1) + rightArea * (numLeafs - i - 1)) / totalArea;
                    if (splitCost < minCost)
                    {
                        minCost    = splitCost;
                        splitIndex = i;
                        splitAxis  = axis;
                        leftBox    = Eigen::AlignedBox<double, TDimensions>(leftBoxes[i].min(), leftBoxes[i].max());
                        rightBox   = Eigen::AlignedBox<double, TDimensions>(rightBoxes[i].min(), rightBoxes[i].max());
                    }
                }
            }

            leftSortedObjects  = new Leaf*[numLeafs];
            rightSortedObjects = new Leaf*[numLeafs];

            leftSortedObjectsList[splitAxis]  = leftSortedObjects;
            rightSortedObjectsList[splitAxis] = rightSortedObjects;
            sortedObjects                     = sortedObjectsList[splitAxis];

            for (i = 0; i <= splitIndex; i++)
            {
                *sortedObjects[i]->isLeft = true;
                leftSortedObjects[i]      = sortedObjects[i];
            }
            for (i = splitIndex + 1; i <= numLeafs - 1; i++)
            {
                *sortedObjects[i]->isLeft              = false;
                rightSortedObjects[i - splitIndex - 1] = sortedObjects[i];
            }

            for (axis = 0; axis < TDimensions; axis++)
            {
                if (axis != splitAxis)
                {
                    leftSortedObjects  = new Leaf*[numLeafs];
                    rightSortedObjects = new Leaf*[numLeafs];

                    leftSortedObjectsList[axis]  = leftSortedObjects;
                    rightSortedObjectsList[axis] = rightSortedObjects;
                    sortedObjects                = sortedObjectsList[axis];

                    leftIndex  = 0;
                    rightIndex = 0;
                    for (i = 0; i < numLeafs; i++)
                    {
                        if (*sortedObjects[i]->isLeft)
                        {
                            leftSortedObjects[leftIndex] = sortedObjects[i];
                            leftIndex++;
                        }
                        else
                        {
                            rightSortedObjects[rightIndex] = sortedObjects[i];
                            rightIndex++;
                        }
                    }
                }
            }
            if (leftIndex == 0)
            {
                printf("BUG left index == 0\n");
                return;
            }
            if (rightIndex == 0)
            {
                printf("BUG right index == 0 left index %d\n", leftIndex);
                return;
            }

            if (leftIndex > 1)
            {
                Node* leftNode = new Node();
                node->addLeft(leftNode);
                createHierarchy(leftNode, leftBox, leftSortedObjectsList, splitIndex + 1);
                node->recomputeBoundingBox();
            }
            else
            {
                node->setLeafs(leftSortedObjects[0], NULL);
                node->recomputeBoundingBox();
            }

            if (rightIndex > 1)
            {
                Node* rightNode = new Node();
                node->addRight(rightNode);
                createHierarchy(rightNode, rightBox, rightSortedObjectsList, numLeafs - splitIndex - 1);
                node->recomputeBoundingBox();
            }
            else
            {
                node->setLeafs(NULL, rightSortedObjects[0]);
                node->recomputeBoundingBox();
            }

            for (int k = 0; k < TDimensions; ++k)
            {
                delete[] leftSortedObjectsList[k];
                delete[] rightSortedObjectsList[k];
            }
            delete[] leftBoxes;
            delete[] rightBoxes;
        }

        /**
         * @brief Creates the bounding boxes on the left.
         * @param leftBoxes Pointer where to store the bounding boxes.
         * @param numShapes Number of shapes.
         * @param sortedObjects Pointer to the sorted shapes.
         */
        void createLeftBoxes(Eigen::AlignedBox<double, TDimensions>* leftBoxes, int numShapes, Leaf** sortedObjects)
        {
            Eigen::AlignedBox<double, TDimensions> box;
            box.setEmpty();
            for (int i = 0; i < numShapes - 1; i++)
            {
                Leaf* leaf = sortedObjects[i];
                box.extend(leaf->shape->worldBounds());
                leftBoxes[i] = Eigen::AlignedBox<double, TDimensions>(box.min(), box.max());
            }
        }

        /**
         * @brief Creates the bounding boxes on the right.
         * @param rightBoxes Pointer where to store the bounding boxes.
         * @param numShapes Number of shapes.
         * @param sortedObjects Pointer to the sorted shapes.
         */
        void createRightBoxes(Eigen::AlignedBox<double, TDimensions>* rightBoxes, int numShapes, Leaf** sortedObjects)
        {
            Eigen::AlignedBox<double, TDimensions> box;
            box.setEmpty();
            for (int i = numShapes - 1; i > 0; i--)
            {
                Leaf* leaf = sortedObjects[i];
                box.extend(leaf->shape->worldBounds());
                rightBoxes[i - 1] = Eigen::AlignedBox<double, TDimensions>(box.min(), box.max());
            }
        }

        /**
         * @brief Computes the circumference of a 2D bounding box.
         * @param bounds Bounding box to compute circumference for.
         * @return Circumference of bounding box.
         */
        double computeArea(const Eigen::AlignedBox2d& bounds) const
        {
            if (bounds.isEmpty())
                return std::numeric_limits<double>::max() / 10;
            Eigen::Vector2d delta = bounds.max() - bounds.min();
            return 2.0 * (delta.x() + delta.y());
        }

        /**
         * @brief Computes the surface area of a 3D bounding box.
         * @param bounds Bounding box to compute surface area for.
         * @return Surface area of bounding box.
         */
        double computeArea(const Eigen::AlignedBox3d& bounds) const
        {
            if (bounds.isEmpty())
                return std::numeric_limits<double>::max() / 10;
            Eigen::Vector3d delta = bounds.max() - bounds.min();
            return 2.0 * (delta.x() * delta.y() + delta.x() * delta.z() + delta.y() * delta.z());
        }

        /**
         * @brief Shapes to add into the tree.
         */
        std::vector<TShape*> mShapes;

        /**
         * @brief Root node of the hierarchy.
         */
        Node* mRoot;

        /**
         * @brief Array of sorted shapes.
         */
        std::array<Leaf**, TDimensions> mSortedObjects;
    };
}
