#define BOOST_TEST_MODULE header-only testPermuteDims
#include <boost/test/included/unit_test.hpp>
#include <mpi.h>
#include "../src/permuteDims.cpp"
#include <cstdint>

int ME = 0, NPROCS = 1;
struct MPISetup
{
    MPISetup()
    {
        int argc = 0;
        char **argv = nullptr;
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &ME);
        MPI_Comm_size(MPI_COMM_WORLD, &NPROCS);
    }
    ~MPISetup() { MPI_Finalize(); }
};

BOOST_TEST_GLOBAL_FIXTURE(MPISetup);


BOOST_AUTO_TEST_CASE(single_proc_small_2d)
{
    std::vector<int> vec(10, 0);
    std::vector<uint64_t> dimLengths = {2, 5};
    std::vector<uint64_t> newDims = {0, 1};

    auto ret = MDPAT::permuteDims<int>(vec, dimLengths, newDims);
    BOOST_TEST(ret == -3);

    newDims = {1, 0};
    ret = MDPAT::permuteDims<int>(vec, dimLengths, newDims);
    BOOST_TEST(ret == 0);

    vec = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    ret = MDPAT::permuteDims<int>(vec, dimLengths, newDims);
    BOOST_TEST(ret == 0);
    BOOST_TEST(vec[0] == 0);
    BOOST_TEST(vec[1] == 5);
    BOOST_TEST(vec[2] == 1);
    BOOST_TEST(vec[3] == 6);
    BOOST_TEST(vec[4] == 2);
    BOOST_TEST(vec[5] == 7);
    BOOST_TEST(vec[6] == 3);
    BOOST_TEST(vec[7] == 8);
    BOOST_TEST(vec[8] == 4);
    BOOST_TEST(vec[9] == 9);
    
}

BOOST_AUTO_TEST_CASE(single_proc_small_3d)
{
    std::vector<int> vec(8, 0);
    std::vector<uint64_t> dimLengths = {2, 2, 2};
    std::vector<uint64_t> newDims = {0, 1, 2};

    auto ret = MDPAT::permuteDims<int>(vec, dimLengths, newDims);
    BOOST_TEST(ret == -3);

    // Start with [ [[0, 1],
    //               [2, 3]],
    //              [[4, 5],    vec[1][0][1] = 5
    //               [6, 7]] ]
    // Get to     [ [[0, 1],
    //               [4, 5]],   vec[0][1][1] = 5
    //              [[2, 3],
    //               [6, 7]] ]
    int val = -1;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 2; ++k)
                vec[i*4+j*2+k] = ++val;
    // vec = {0, 1, 2, 3, 4, 5, 6, 7};
    newDims[0] = 1;
    newDims[1] = 0;
    
    ret = MDPAT::permuteDims<int>(vec, dimLengths, newDims);
    BOOST_TEST(ret == 0);

    val == -1;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 2; ++k)
                BOOST_TEST(vec[j*4+i*2+k] = ++val);
}

BOOST_AUTO_TEST_CASE(two_proc_2d)
{
    std::vector<int> vec(8, 0);
    std::vector<unsigned long> dimLengths = {2, 4};
    std::vector<unsigned long> newDims = {1, 0};

    int val = -1;
    for (int i = 0; i < 8; ++i)
        vec[i] = ++val + ME*8;

    // Before:
    // Proc 0: [ [ 0,  1,  2,  3],
    //           [ 4,  5,  6,  7] ]
    // Proc 1: [ [ 8,  9, 10, 11],
    //           [12, 13, 14, 15] ]

    // After:
    // Proc 0: [ [ 0,  4,  8, 12],
    //           [ 1,  5,  9, 13] ]
    // Proc 1: [ [ 2,  6, 10, 14],
    //           [ 3,  7, 11, 15] ]

    auto ret = MDPAT::permuteDimsParallel<int>(vec, dimLengths, newDims, ME, NPROCS, MPI_COMM_WORLD);
    BOOST_TEST(ret == 0);

    if (ME == 0)
    {
        BOOST_TEST(vec[0] == 0);
        BOOST_TEST(vec[1] == 4);
        BOOST_TEST(vec[2] == 8);
        BOOST_TEST(vec[3] == 12);
        BOOST_TEST(vec[4] == 1);
        BOOST_TEST(vec[5] == 5);
        BOOST_TEST(vec[6] == 9);
        BOOST_TEST(vec[7] == 13);
    }
    else
    {
        BOOST_TEST(vec[0] == 2);
        BOOST_TEST(vec[1] == 6);
        BOOST_TEST(vec[2] == 10);
        BOOST_TEST(vec[3] == 14);
        BOOST_TEST(vec[4] == 3);
        BOOST_TEST(vec[5] == 7);
        BOOST_TEST(vec[6] == 11);
        BOOST_TEST(vec[7] == 15);
    }
}
