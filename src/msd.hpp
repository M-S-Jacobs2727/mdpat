#include <cstdint>
#include <filesystem>
#include <vector>

#include <mpi.h>
#include <omp.h>

#include "permuteDims.hpp"
#include "readTrajectories.hpp"
#include "selectFromArray.hpp"
#include "splitValues.hpp"

namespace fs = std::filesystem;

namespace MDPAT
{
    template <typename T>
    struct msd_results
    {
        std::vector<uint64_t> timegaps;
        std::vector<T> msd;
    };

    template <typename T>
    msd_results<T> meanSquaredDisplacement(
        std::vector<int32_t> &typelist,
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

    void meanSquaredDisplacement(
        fs::path outfile,
        fs::path directory,
        StepRange stepRange,
        double timestep,
        uint32_t atomType,
        uint64_t minGap, // in number of steps
        uint64_t maxGap, // in number of steps
        uint32_t dim,    // number of spatial dimensions
        int me,
        int nprocs,
        MPI_Comm comm);

    void meanSquaredDisplacementOMP(
        fs::path outfile,
        fs::path directory,
        StepRange stepRange,
        double timestep,
        uint32_t atomType,
        uint64_t minGap, // in number of steps
        uint64_t maxGap, // in number of steps
        uint32_t dim,    // number of spatial dimensions
        int me,
        int nprocs,
        MPI_Comm comm);

}
