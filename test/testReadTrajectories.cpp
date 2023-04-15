#define BOOST_TEST_MODULE testReadTrajectories
#include <boost/test/unit_test.hpp> 
#include "readTrajectories.hpp"

BOOST_AUTO_TEST_CASE(single_text_file) {
    const std::string filename = "test/testFiles/dump.0";
    std::ifstream dumpFile(filename);
    BOOST_TEST(dumpFile.good());
    
    std::vector<int> columnFlag;
    columnFlag.resize(4, 1);
    
    std::vector<Frame<float>*> frames;
    frames.resize(1);
    frames[0] = new Frame<float>;
    
    auto err_code = readDumpTextFile<float>(dumpFile, columnFlag, frames);

    BOOST_TEST(err_code == 0);
    
    BOOST_TEST(frames[0]->nAtoms == 500L);
    BOOST_TEST(frames[0]->timestep == 0L);

    BOOST_TEST(frames[0]->box[0] == 0.0F);
    BOOST_TEST(frames[0]->box[1] == 8.3979809569125372F);
    BOOST_TEST(frames[0]->box[2] == 0.0F);
    BOOST_TEST(frames[0]->box[3] == 8.3979809569125372F);
    BOOST_TEST(frames[0]->box[4] == 0.0F);
    BOOST_TEST(frames[0]->box[5] == 8.3979809569125372F);

    BOOST_TEST(frames[0]->atoms.size() == 2000);

    // Atom 1
    BOOST_TEST(frames[0]->atoms[0] == 1.0F);
    BOOST_TEST(frames[0]->atoms[1] == 0.0F);
    BOOST_TEST(frames[0]->atoms[2] == 0.0F);
    BOOST_TEST(frames[0]->atoms[3] == 0.0F);

    // Atom 6 (the first out-of-order atom in the file)
    BOOST_TEST(frames[0]->atoms[20] == 1.0F);
    BOOST_TEST(frames[0]->atoms[21] == 0.3F);
    BOOST_TEST(frames[0]->atoms[22] == 0.1F);
    BOOST_TEST(frames[0]->atoms[23] == 0.0F);
}
