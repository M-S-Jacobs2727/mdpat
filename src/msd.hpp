#include <cstdint>
#include <filesystem>
#include <vector>

#include "mpi.h"

#include "readTrajectories.hpp"
#include "permuteDims.hpp"
#include "selectFromArray.hpp"
#include "splitValues.hpp"

namespace fs = std::filesystem;

namespace MDPAT
{
template<typename T>
struct msd_results{
    std::vector<uint64_t> timegaps;
    std::vector<T> msd;
};

template<typename T>
msd_results<T> msd(std::vector<int32_t> typelist,
                   std::vector<T> unwrapped_trajectories,
                   const uint64_t first_frame,
                   const uint64_t num_frames,
                   const uint64_t num_atoms,
                   const uint64_t num_spatial_dims,
                   const uint64_t gap_start,
                   const uint64_t gap_end,
                   const int me,
                   const int nprocs,
                   MPI_Comm comm);

void msd(fs::path outfile,
         fs::path directory,
         StepRange stepRange,
         double timestep,
         uint32_t atomType,
         uint64_t minGap,
         uint64_t maxGap,
         uint32_t dim);

}
