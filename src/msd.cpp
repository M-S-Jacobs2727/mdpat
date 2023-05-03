#include "msd.hpp"

namespace MDPAT
{
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
        MPI_Comm comm)
    {
        // * As a whole, this library should read the trajectories once and send
        // * them whole-cloth to each method (such as msd). However, since a method
        // * such as msd will alter the trajectories, and copying may take up too
        // ? much memory, should I instead just read the trajectories a few times?

        // Expects unwrapped_trajectories to come in as a 1D vector representing 3D
        // data of the shape (num_frames, num_atoms, num_spatial_dims)

        // Using permuteDims, it switches the order to (num_atoms, num_spatial_dims, num_frames),
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

        std::vector<uint64_t> dimLengths = {num_frames, num_atoms, num_spatial_dims};
        const std::vector<uint64_t> newDims = {2, 0, 1};

        // TODO: select atoms from typelist

        auto ret = permuteDimsParallel<T>(unwrapped_trajectories, dimLengths, newDims, me, nprocs, comm);
        if (ret)
            throw -1;

        auto my_natoms = dimLengths[0];
        auto total_num_frames = dimLengths[2];

        msd_results<T> results;
        uint64_t results_size = gap_end - gap_start + 1;
        results.timegaps.resize(results_size);
        results.msd.resize(results_size);
        T tmp;

        for (uint64_t gap = gap_start; gap <= gap_end; ++gap)
            results.timegaps[gap - gap_start] = gap;

        for (uint64_t atom = 0; atom < my_natoms; ++atom)
        {
            for (uint64_t col = 0; col < num_spatial_dims; ++col)
            {
                for (uint64_t gap = gap_start; gap <= gap_end; ++gap)
                {
                    for (uint64_t frame = 0; frame + gap < total_num_frames; ++frame)
                    {
                        tmp = unwrapped_trajectories[atom * num_spatial_dims * total_num_frames + col * total_num_frames + frame + gap] -
                              unwrapped_trajectories[atom * num_spatial_dims * total_num_frames + col * total_num_frames + frame];
                        results.msd[gap - gap_start] += tmp * tmp;
                    }
                }
            }
        }

        return results;
    }

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
        MPI_Comm comm)
    {
        uint64_t nSteps = (stepRange.endStep - stepRange.initStep) / stepRange.dumpStep + 1;
        auto [myFirstStep, myNumSteps] = splitValues(nSteps, me, nprocs);

        StepRange myStepRange(myFirstStep,
                              myFirstStep + stepRange.dumpStep * (myNumSteps - 1),
                              stepRange.dumpStep);

        auto atoms = getTrajectories(false, myStepRange, directory, atomType, dim);
        auto nAtoms = atoms.size() / dim / myNumSteps;

        std::vector<uint64_t> dimLengths = {myNumSteps, nAtoms, dim};
        std::vector<uint64_t> newDims = {1, 2, 0};
        permuteDimsParallel(atoms, dimLengths, newDims, me, nprocs, comm);
        auto myNumAtoms = atoms.size() / dim / nSteps;

        minGap /= myStepRange.dumpStep;
        maxGap /= myStepRange.dumpStep;
        uint64_t numGaps = maxGap - minGap + 1;

        double rsq, dx;
        std::vector<double> msd(numGaps, 0.0);

        for (uint64_t atom = 0; atom < myNumAtoms; ++atom)
        {
            for (uint64_t col = 0; col < dim; ++col)
            {
                for (uint64_t gap = minGap; gap <= maxGap; ++gap)
                {
                    rsq = 0;
                    for (uint64_t frame = 0; frame + gap < nSteps; ++frame)
                    {
                        dx = atoms[atom * dim * nSteps + col * nSteps + frame + gap] -
                             atoms[atom * dim * nSteps + col * nSteps + frame];
                        rsq += dx * dx;
                    }
                    msd[gap - minGap] += rsq;
                }
            }
        }

        if (me)
            MPI_Reduce(msd.data(), msd.data(), msd.size(), MPI_DOUBLE, MPI_SUM, 0, comm);
        else
            MPI_Reduce(MPI_IN_PLACE, msd.data(), msd.size(), MPI_DOUBLE, MPI_SUM, 0, comm);

        if (me == 0)
        {
            std::ofstream outstream(outfile);
            for (uint64_t gap = minGap; gap <= maxGap; ++gap)
                outstream << gap * myStepRange.dumpStep * timestep << ' '
                          << msd[gap - minGap] / nAtoms / (nSteps - gap) << '\n';
            outstream.close();
        }
    }

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
        MPI_Comm comm)
    {
        uint64_t nSteps = (stepRange.endStep - stepRange.initStep) / stepRange.dumpStep + 1;
        auto [myFirstStep, myNumSteps] = splitValues(nSteps, me, nprocs);

        StepRange myStepRange(myFirstStep,
                              myFirstStep + stepRange.dumpStep * (myNumSteps - 1),
                              stepRange.dumpStep);

        auto atoms = getTrajectories(false, myStepRange, directory, atomType, dim);
        auto nAtoms = atoms.size() / dim / myNumSteps;

        std::vector<uint64_t> dimLengths = {myNumSteps, nAtoms, dim};
        std::vector<uint64_t> newDims = {1, 2, 0};
        permuteDimsParallel(atoms, dimLengths, newDims, me, nprocs, comm);
        auto myNumAtoms = atoms.size() / dim / nSteps;

        minGap /= myStepRange.dumpStep;
        maxGap /= myStepRange.dumpStep;
        uint64_t numGaps = maxGap - minGap + 1;

        double rsq, dx;
        std::vector<double> msd(numGaps, 0.0);

#pragma omp distribute collapse(3)
        for (uint64_t atom = 0; atom < myNumAtoms; ++atom)
        {
            for (uint64_t col = 0; col < dim; ++col)
            {
                for (uint64_t gap = minGap; gap <= maxGap; ++gap)
                {
                    rsq = 0;
#pragma omp simd reduction(+ : rsq)
                    for (uint64_t frame = 0; frame < nSteps - gap; ++frame)
                    {
                        dx = atoms[atom * dim * nSteps + col * nSteps + frame + gap] -
                             atoms[atom * dim * nSteps + col * nSteps + frame];
                        rsq += dx * dx;
                    }
#pragma omp atomic
                    msd[gap - minGap] += rsq;
                }
            }
        }

        if (me)
            MPI_Reduce(msd.data(), msd.data(), msd.size(), MPI_DOUBLE, MPI_SUM, 0, comm);
        else
            MPI_Reduce(MPI_IN_PLACE, msd.data(), msd.size(), MPI_DOUBLE, MPI_SUM, 0, comm);

        if (me == 0)
        {
            std::ofstream outstream(outfile);
            for (uint64_t gap = minGap; gap <= maxGap; ++gap)
                outstream << gap * myStepRange.dumpStep * timestep << ' '
                          << msd[gap - minGap] / nAtoms / (nSteps - gap) << '\n';
            outstream.close();
        }
    }

}
