#define BOOST_TEST_MODULE header-only testPermuteDims
#include <boost/test/included/unit_test.hpp>
#include <filesystem>
#include "../src/readTrajectories.hpp"
#include "../src/msd.hpp"

namespace fs = std::filesystem;

BOOST_AUTO_TEST_CASE(msd_serial_1)
{
    fs::path outfile("msd_test.txt");
    fs::path directory("../test/testFiles");
    MDPAT::StepRange stepRange(0UL, 250UL, 50UL);
    double timestep = 0.005;
    uint32_t atomType = 1;
    uint64_t minGap = 0UL;
    uint64_t maxGap = 100UL;
    uint32_t dim = 3;
    int me = 0;
    int nprocs = 1;
    
    MDPAT::meanSquaredDisplacement(
        outfile,
        directory,
        stepRange,
        timestep,
        atomType,
        minGap, // in number of steps
        maxGap, // in number of steps
        dim,    // number of spatial dimensions
        me,
        nprocs,
        MPI_COMM_WORLD);
}
