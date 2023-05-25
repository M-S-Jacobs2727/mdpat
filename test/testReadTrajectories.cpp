#define BOOST_TEST_MODULE header-only testReadTrajectories
#include <boost/test/included/unit_test.hpp> 
#include "../src/readTrajectories.hpp"
#include <filesystem>

namespace fs = std::filesystem;

BOOST_AUTO_TEST_CASE(find_coord_column)
{
    MDPAT::Frame <float>frame;
    frame.columnLabels = {"id", "type", "xs", "ys", "zs"};

    auto ret = MDPAT::findCoordColumn(frame, 'x', true);
    BOOST_TEST(ret == -2);

    ret = MDPAT::findCoordColumn(frame, 'y', true);
    BOOST_TEST(ret == -3);

    ret = MDPAT::findCoordColumn(frame, 'z', true);
    BOOST_TEST(ret == -4);
}

BOOST_AUTO_TEST_CASE(get_trajectories)
{
    MDPAT::StepRange steprange(0UL, 250UL, 50UL);

    BOOST_TEST(steprange.nSteps == 6);

    fs::path directory("../test/testFiles");
    auto atoms = MDPAT::getTrajectories(directory, steprange, true, 0, 3);

    BOOST_TEST(atoms.size() == 500*6*3  );

    BOOST_TEST(atoms[0] == 0.0F);
    BOOST_TEST(atoms[1] == 0.0F);
    BOOST_TEST(atoms[2] == 0.0F);
    BOOST_TEST(atoms[3] == 0.83979809569125372F);
    BOOST_TEST(atoms[4] == 0.83979809569125372F);
    BOOST_TEST(atoms[5] == 0.0F);
}
