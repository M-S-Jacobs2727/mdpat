cmake_minimum_required(VERSION 3.10)
project(MDPAT_TEST VERSION 0.1.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR})

find_package(Boost COMPONENTS system filesystem unit_test_framework REQUIRED)
add_executable(test
    ${TEST_DIR}/testReadTrajectories.cpp
)
target_link_libraries(test
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)