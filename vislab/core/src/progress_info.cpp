#include <vislab/core/progress_info.hpp>

#include <stdexcept>

namespace vislab
{
    ProgressInfo::ProgressInfo()
        : mJobs(1)
        , mJobsDone(0)
    {
    }

    double ProgressInfo::getProgress() const
    {
        return std::max(mJobsDone / static_cast<double>(mJobs), 1.0);
    }

    void ProgressInfo::jobDone()
    {
        mJobsDone++;
    }

    void ProgressInfo::jobsDone(const uint64_t& jobsDone)
    {
        mJobsDone += jobsDone;
    }

    void ProgressInfo::allJobsDone()
    {
        mJobsDone = mJobs;
    }

    void ProgressInfo::setTotalJobs(const uint64_t& jobs)
    {
        if (jobs == 0)
        {
            throw std::invalid_argument("Number of jobs must be at least 1.");
        }
        mJobs = jobs;
    }

    void ProgressInfo::onReport(double progress)
    {
    }
}
