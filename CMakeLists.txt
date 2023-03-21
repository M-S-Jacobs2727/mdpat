cmake_minimum_required(VERSION 3.10)
project(MDPAT VERSION 0.1.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(MDPAT_DIR ${CMAKE_SOURCE_DIR})
set(MDPAT_SRC_DIR ${MDPAT_DIR}/src)
set(MDPAT_TOOLS_DIR ${MDPAT_DIR}/tools)
set(MDPAT_TEST_DIR ${MDPAT_DIR}/test)

add_executable(mdpat
    ${MDPAT_SRC_DIR}/main.cpp
    ${MDPAT_SRC_DIR}/bcastContainers.cpp
    ${MDPAT_SRC_DIR}/initNode.cpp
    ${MDPAT_SRC_DIR}/permuteDims.cpp
    ${MDPAT_SRC_DIR}/readInput.cpp
    ${MDPAT_SRC_DIR}/readTrajectories.cpp
    ${MDPAT_SRC_DIR}/splitValues.cpp
    ${MDPAT_SRC_DIR}/selectFromArray.cpp
)

add_executable()
