#include <array>
#include <vector>
#include <cstdint>

#include "mpi.h"

#include "permuteDims.hpp"

template<typename T>
struct msd_results{
    std::vector<uint32_t> timegaps;
    std::vector<T> msd;
};

template<typename T>
msd_results<T> msd(std::vector<int32_t> typelist,
                   std::vector<T> unwrapped_trajectories,
                   uint32_t first_frame,
                   uint32_t num_frames,
                   uint32_t num_atoms,
                   uint32_t num_cols,
                   uint32_t gap_start,
                   uint32_t gap_end)
{
    // Expects unwrapped_trajectories to come in as a 1D vector representing 3D
    // data of the shape (num_frames, num_atoms, num_cols)

    // Using permuteDims, it switches the order to (num_atoms, num_cols, num_frames),
    // because the MSD calculation can take place over each atom and spatial dimension
    // independently. Thus, the GPU can compute the MSD value for one atom and spatial
    // dimension per thread and accumulate the results.

    // The issue is that this function is called by each process with its own set of
    // frames. So after permuting the dimensions, we then have to reorganize the data
    // among the processes to give each proc all the frames corresponding to its atoms
    // and spatial dimensions.

    // Changing my mind: I'll do what I did before and send the data to proc 0 to 
    // reorganize and redistribute. It's the easiest, and if I want to make it better
    // later on, I can change it.

    MPI_Gatherv(&unwrapped_trajectories[0], unwrapped_trajectories.size(), MPI_FLOAT, );


    const std::array<uint32_t, 3> dimLengths = {num_frames, num_atoms, num_cols};
    const std::array<uint32_t, 3> newDims = {2, 3, 1};
    auto new_trajectories = permuteDims(unwrapped_trajectories, dimLengths, newDims);
}
