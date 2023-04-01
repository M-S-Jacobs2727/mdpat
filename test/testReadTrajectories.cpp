#define BOOST_TEST_MODULE testReadTrajectories
#include <boost/test/unit_test.hpp> 
#include "readTrajectories.hpp"

BOOST_AUTO_TEST_CASE(single_text_file) {
    Dump<float> dump;
    std::vector<int> columnFlag;
    columnFlag.resize(4, 1);
    const std::string dumpfile = "test/testFiles/dump.0";
    auto err_code = readDumpText<float>(dumpfile, 500, columnFlag, dump);
    BOOST_TEST(err_code == 0);
}
