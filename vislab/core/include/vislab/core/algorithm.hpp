#pragma once

#include "data.hpp"
#include "input_port.hpp"
#include "output_port.hpp"
#include "parameter.hpp"
#include "progress_info.hpp"
#include "update_info.hpp"

namespace vislab
{
    /**
     * @brief Base class for algorithms.
     */
    class Algorithm : public Object
    {
        VISLAB_INTERFACE(Algorithm, Object)

    public:
        /**
         * @brief Constructor.
         */
        Algorithm();

        /**
         * @brief Destructor.
         */
        virtual ~Algorithm();

        /**
         * @brief Invokes the execution of the algorithm.
         * @param progress Can be used to monitor the progress.
         * @return Information object with potential error message.
         */
        UpdateInfo update(ProgressInfo& progress);

        /**
         * @brief Invokes the execution of the algorithm.
         * @return Information object with potential error message.
         */
        UpdateInfo update();

        /**
         * @brief Flag that determines whether this algorithm is currently active.
         */
        BoolParameter active;

    protected:
        /**
         * @brief Internal computation function
         * @param progress Optional progress info.
         * @return Information about the completion of the computation, including a potential error message.
         */
        virtual UpdateInfo internalUpdate(ProgressInfo& progress) = 0;
    };

#define VISLAB_ALGORITHM(CLASS, ...)                                                          \
    VISLAB_OBJECT(CLASS, __VA_ARGS__)
}
