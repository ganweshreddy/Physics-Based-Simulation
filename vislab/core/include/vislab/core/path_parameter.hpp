#pragma once

#include "parameter.hpp"

namespace vislab
{
    /**
     * @brief Parameter that represents the path to a file.
     */
    class PathParameter : public Parameter<std::string>
    {
        VISLAB_OBJECT(PathParameter, Parameter<std::string>)

    public:
        /**
         * @brief Flag that determines whether the path is to be read or written to.
         */
        enum class EFile
        {
            /**
             * @brief Read from the path.
             */
            In,

            /**
             * @brief Write to the path.
             */
            Out
        };

        /**
         * @brief Constructor.
         */
        PathParameter();

        /**
         * @brief Constructor.
         */
        PathParameter(const EFile& fileDirection, const std::string& filter);

        /**
         * @brief Gets the file reading direction.
         * @return File reading direction.
         */
        const EFile& getFileDirection() const;

        /**
         * @brief Sets the file reading direction.
         * @param file New file reading direction to set.
         */
        void setFileDirection(const EFile& file);

        /**
         * @brief Gets a filter string for the UI, which restricts to certain files with a certain extension.
         * @return Filter string for the UI.
         */
        const std::string& getFilter() const;

        /**
         * @brief Sets a filter string for the UI, which restricts to certain files with a certain extension.
         * @param filter New filter string to set.
         */
        void setFilter(const std::string& filter);

        /**
         * @brief Tests if another parameter stores the exact same values and whether they have the same visibilty.
         * @param other Other parameter to compare with.
         * @return True, if the values and visibility are the same.
         */
        bool operator==(const PathParameter& other) const;

    private:
        /**
         * @brief Filter string, which controls the file extensions that can be selected.
         */
        std::string mFilter;

        /**
         * @brief File reading direction, which is in or out.
         */
        EFile mFile;
    };
}
